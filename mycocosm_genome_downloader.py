#!/usr/bin/env python

"""
This script downloads all fungal genomes from JGI: assembly (fasta) + annotation (gff)

Also, get a summary of metadata:
TaxId   JGI_short_name  Complete name   Original name   +files location
(complete name must drop the "v1.0" etc.)
Note that JGI only lists the latest version (i.e. if v2.0 is available for a 
certain assembly, v1.0 will not be available anymore)
    
    
generic XML quick guide (to manipulate with etree):
This is a node:
    <data year="2018">
        1234
    </data>

'data' = node.tag
'year' = node.attrib (dictionary)
'1234' = node.text

Elements can be nested.

Version 2:
requires: lxml etetoolkit ete3 ete_toolchain

#TODO: exit if JGI logging in is not successful
"""

import sys
import os
import argparse
from pathlib import Path
import csv
import subprocess as sp
from subprocess import DEVNULL, STDOUT, CalledProcessError
from lxml import etree
from getpass import getpass
from ete3 import NCBITaxa
import codecs # there's some strange character from MycoCosm's list
import time
from datetime import datetime
from collections import defaultdict
import shutil

__author__ = "Jorge Navarro"
__email__ = "j.navarro@wi.knaw.nl"
__version__ = "v1.0"


class JGI_Project:
    """
    This class keeps information about each JGI sequencing project
    """
    
    def __init__(self):
        self.TaxId = ""             # NCBI TaxId (not lineage)
        self.lineage_set = set()    # to quickly check major fungal taxonomy branches
        self.lineage_list = []

        self.portal = ""            # Fam[:2] + Sp[:1] + mystery number. a.k.a.: "shortname"
        self.name = ""              # Portal full name. Include things like "v2.0"
        self.org_name = ""          # Portal name without strings like "v2.0"
        
        self.assembly_file = ""
        self.assembly_url = ""
        self.assembly_size = 0
        self.gff_file = ""
        self.gff_url = ""
        self.gff_timestamp = ""
        self.gff_size = 0
        self.proteome_file = ""
        self.proteome_url = ""

        self.project_path = ""      # Path object

        self.date = ""              # Ideally, a Time object
        
        return


def command_parser():
    parser = argparse.ArgumentParser(description=f"Downloads fungal genomes from JGI's MycoCosm. Version: {__version__}")

    group_get = parser.add_argument_group("Get necessary data")
    group_get.add_argument("-u", "--update", help="Updates NCBI Taxonomy database\
        and exits. Used to infer lineage from TaxId", default=False, 
        action="store_true")
    group_get.add_argument("--getgenomelist", help="Gets a fresh download \
        of MycoCosm's genome list 'MycoCosm_Genome_list.csv' and exits. Used to\
        get metadata of each project (TaxId, full name). Stored in the folder \
        specified by --outputfolder", default=False, action="store_true")
    group_get.add_argument("--getxml", help="Gets a fresh download of JGI's xml\
        file and exits. Use to obtain routes to each file. Stored in the folder \
        specified by --outputfolder", default=False, action="store_true")
    group_get.add_argument("--getprevious", type=Path, help="Gets a list of\
        locations of .gz files, downloaded from a previous run. Argument to \
        this parameter is the base folder to start looking recursively for \
        .gz files. File will be stored in the same place as this script (name:\
        'previously_downloaded_files.tsv'")

    group_input = parser.add_argument_group("Input")
    group_input.add_argument("-c", "--csv", help="List of genomes in JGI with \
        some metadata like TaxId", type=Path)
    group_input.add_argument("-x", "--xml", help="JGI file with links to files",
        type=Path)
    group_input.add_argument("-p", "--previous", help="File listing the \
        locations of previously-downloaded files to copy from, instead of \
        re-downloading", type=Path)

    group_output = parser.add_argument_group("Output")
    group_output.add_argument("-o", "--outputfolder", help="Base folder to  \
        download all files", type=Path, 
        default=(Path(__file__).parent / "output"))
    group_output.add_argument("-s", "--simulate", help="Don't download \
        anything, just create the folder structure and taxonomy file", \
        default=False, action="store_true")

    return parser.parse_args()




#*****************************************************************************

def JGI_login(cookie_path):
    """
    Get (physical) cookie. Needed to download files from JGI
    """

    print("Please enter your JGI credentials:")
    user = input("Username: ")
    password = getpass()

    if user == "" or password == "":
        sys.exit("Missing info for JGI login")
    
    # This should be equivalent to
    # curl 'https://signon.jgi.doe.gov/signon/create' --data-urlencode 'login=USER_NAME' --data-urlencode 'password=USER_PASSWORD' -c cookies > /dev/null
    # See: https://genome.jgi.doe.gov/portal/help/download.jsf#api
    command = []
    command.append("curl")
    command.append("https://signon.jgi.doe.gov/signon/create")
    command.append("--data-urlencode")
    command.append(f"login={user}")
    command.append("--data-urlencode")
    command.append(f"password={password}")
    command.append("-c")
    command.append(f"{cookie_path}")

    print(" Logging in:")
    print("--------------------------------------------")
    proc = sp.run(command, shell=False, stdout=DEVNULL)
    try:
        proc.check_returncode()
    except sp.CalledProcessError:
        print("Error logging in: {}".format(proc.stderr))
        return False
    else:
        print("--------------------------------------------\n")
        return True


def get_JGI_genome_list(outputfolder):
    cmd = []
    cmd.append("curl")
    cmd.append("https://mycocosm.jgi.doe.gov/ext-api/mycocosm/catalog/download-group?flt=&seq=all&pub=all&grp=fungi&srt=released&ord=asc")
    cmd.append("-o")
    cmd.append(str(outputfolder/"MycoCosm_Genome_list.csv"))

    print(" Downloading genome list")
    print("--------------------------------------------")
    proc = sp.run(cmd, stderr=STDOUT, encoding="utf-8")
    try:
        proc.check_returncode()
    except sp.CalledProcessError:
        print("Error: {}".format(proc.stderr))
        return False
    else:
        print("--------------------------------------------\n")
        return True


def remove_version(label):
    '''
    Try to remove assembly version from label
    '''
    suffix = label.strip().split(" ")[-1]
    if (suffix[0] == "v" or suffix[0] == "V") and suffix[1].isdigit() and suffix[-1].isdigit():
        # Assume we have one of 
        # v3, v1.0, v2.0, v4.0, v1.2, v2.2, v1.1, v3.0, v2, etc.
        return " ".join(label.strip().split(" ")[:-1])
    else:
        return label


# main branches of the JGI tree. See https://genome.jgi.doe.gov/programs/fungi/index.jsf
jgi_tree_branches = set({"pucciniomycotina", 
                    "ustilaginomycotina", 
                    "agaricomycetes", 
                    "dacrymycetes", 
                    "tremellomycetes", 
                    "wallemiomycetes", 
                    "pezizomycetes", 
                    "orbiliomycetes", 
                    "eurotiomycetes", 
                    "dothideomycetes", 
                    "lecanoromycetes", 
                    "leotiomycetes", 
                    "sordariomycetes", 
                    "xylonomycetes", 
                    "saccharomycotina", 
                    "taphrinomycotina", 
                    "glomeromycotina", 
                    "mortierellomycotina", 
                    "mucoromycotina", 
                    "zoopagomycotina", 
                    "entomophthoromycotina", 
                    "kickxellomycotina", 
                    "blastocladiomycota", 
                    "chytridiomycetes", 
                    "monoblepharidomycetes", 
                    "neocallimastigomycota", 
                    "microsporidia", 
                    "cryptomycota"})
def get_final_output_folder(o, lineage_set):
    """
    Assembles the path that will contain the files for a particular organism
    based on its lineage
    """
    
    # build folder structure
    output_folder = o
    if "dikarya" in lineage_set:
        output_folder = output_folder / "DIKARYA"
    
        if "ascomycota" in lineage_set:
            output_folder = output_folder / "ASCOMYCOTA"
            
            if "pezizomycotina" in lineage_set:
                output_folder = output_folder / "PEZIZOMYCOTINA"
                
        elif "basidiomycota" in lineage_set:
            output_folder = output_folder / "BASIDIOMYCOTA"
            
            if "agaricomycotina" in lineage_set:
                output_folder = output_folder / "AGARICOMYCOTINA"
                
    else:
        output_folder = output_folder / "no_rank"
        
        if "mucoromycota" in lineage_set:
            output_folder = output_folder / "MUCOROMYCOTA"
            
        if "zoopagomycota" in lineage_set:
            output_folder = output_folder / "ZOOPAGOMYCOTA"
        
        if "chytridiomycota" in lineage_set:
            output_folder = output_folder / "CHYTRIDIOMYCOTA" / "no_rank"
        
    # build last folder
    if len(lineage_set & jgi_tree_branches) == 1:
        output_folder = output_folder / list(lineage_set & jgi_tree_branches)[0].upper()
    else:
        output_folder = output_folder / "no_rank"
        
    return output_folder


# Translates project shortname to organism name. Some had issues with 
# the encoding
manual_project_names = {"Boledp1":"Boletus edulis Přilba v1.0",
                        "Rusoch1":"Russula ochroleuca Přilba v1.0",
                        "Amarub1":"Amanita rubescens Přilba v1.0",
                        "Ruseme1":"Russula emetica Přilba v1.0"}
# metaprojects or old versions we don't want to keep
portals_to_remove = {"Rhoto_IFO0880_2", 
                     "Aciri1_meta",
                     "Pospl1" # superseded by PosplRSB12_1
                     }
def read_mycocosm_csv(csvpath):
    orgdict = dict()
    # MycoCosm file is downloaded with ISO-8859-1 encoding. Use this to avoid UnicodeDecodeError, see 
    # https://stackoverflow.com/questions/12468179/unicodedecodeerror-utf8-codec-cant-decode-byte-0x9c
    with codecs.open(csvpath, "r", encoding="utf-8", errors="replace") as csvfile:
    #with open(csvpath, "r", encoding="utf-8", errors="replace") as csvfile:
        # csv reader is needed because the "name" column contains items that 
        # contain commas, so a simple split() is not enough
        csvreader = csv.DictReader(csvfile)
        for l in csvreader:
            fungus = JGI_Project()
            taxid = l['NCBI Taxon']
            shortname = l['portal']
            fullname = l['name']

            # These taxids don't work for some reason and updating NCBI's
            # database doesn't fix it, so they have to be fixed manually
            # mostly because of "this entry was merged into another entry"
            alt_taxid = {"1140396":"2735509",
                        "585595":"2587404",
                        "5145":"2587412",
                        "1437435":"2211651"}

            if taxid in alt_taxid:
                taxid = alt_taxid[taxid]

            fungus.portal = shortname
            fungus.TaxId = taxid
            if shortname in manual_project_names:
                fungus.name = manual_project_names[shortname]
            else:
                fungus.name = fullname
            fungus.org_name = remove_version(fullname)

            ncbi = NCBITaxa()
            try:
                ncbi_tax_lineage_ids = ncbi.get_lineage(taxid)
            except ValueError:
                ncbi_tax_lineage_ids = []
                print("Can't find TaxId for {} ({}) with ete (try using --update)".format(taxid, shortname))
            
            ncbi_tax_names = ncbi.get_taxid_translator(ncbi_tax_lineage_ids)
            lineage_set = set( [ncbi_tax_names[tid].lower() for tid in ncbi_tax_lineage_ids] )
            lineage_list_str = ",".join([ncbi_tax_names[tid].lower() for tid in ncbi_tax_lineage_ids])

            fungus.lineage_set = lineage_set
            fungus.lineage_list = lineage_list_str

            fungus.project_path = get_final_output_folder(Path("./"), lineage_set)

            orgdict[shortname] = fungus

    # Remove projects we don't want to include
    for d in portals_to_remove:
        if d in orgdict:
            del orgdict[d]

    return orgdict


def get_JGI_xml(outputfolder, cookie_path):
    # Get list of files. Should be equivalent to
    #curl 'http://genome.jgi.doe.gov/portal/ext-api/downloads/get-directory?organism=PhytozomeV10' -b cookies > get-directory
    
    if not cookie_path.is_file():
        sys.exit("Error (get_JGI_xml): cookie not a valid file")
    
    target="fungi" # change target here to make tests
    # target = "Gyrinf1"

    command = []
    command.append("curl")
    command.append("https://genome.jgi.doe.gov/portal/ext-api/downloads/get-directory?organism={}".format(target))
    command.append("-b")
    command.append(str(cookie_path))
    
    print(" Downloading list of files: ")
    # Can't seem to be able to capture the output of the sp.run() command and
    # send it directly pretty-printed to the file; using an intermediate file
    xmlfile = o / "MycoCosm_data.xml"
    # Note: due to the requirements of the environmnent's packages, we're
    # forced to use python 3.6. More recent versions have more options for 
    # subprocess.run(), see https://docs.python.org/3.6/library/subprocess.html
    with open(xmlfile, "w") as f:
        print(" Downloading xml file list")
        print("--------------------------------------------")
        sp.run(command, encoding="utf-8", stdout=f)
        print("--------------------------------------------\n")
    
    print(" Re-formatting xml file")
    parser = etree.XMLParser(remove_blank_text=True)
    xml = etree.parse(str(xmlfile), parser)
    with open(xmlfile, "wb") as f:
        f.write(etree.tostring(xml, pretty_print=True))
    
    return True
    

# Ignore these filenames as they're not related to assemblies, are old versions or
# are assemblies of meta-samples
exclude_assemblies = set(["1034997.Tuber_borchii_Tbo3840.standard.main.scaffolds.fasta.gz",
                         "Spofi1.draft.mito.scaffolds.fasta.gz", 
                         "Patat1.draft.mito.scaffolds.fasta.gz",
                         "PleosPC9_1_Assembly_scaffolds.fasta.gz",
                         "Neuhi1_PlasmidAssemblyScaffolds.fasta.gz",
                         "CocheC5_1_assembly_scaffolds.fasta.gz",
                         "Alternaria_brassicicola_masked_assembly.fasta.gz",
                         "Aciri1_meta_AssemblyScaffolds.fasta.gz",
                         "Rhoto_IFO0880_2_AssemblyScaffolds.fasta.gz",
                         "StenotrophomonasSp_AssemblyScaffolds.fasta.gz",
                         "PseudomonasSp_AssemblyScaffolds.fasta.gz"
                         ])
def annotate_assembly(organisms_csv, xml_Assembly):    
    duplicated_names = defaultdict(list)
    # NOTE:
    # projects' names differ between what's annotated in the csv and xml files
    # for example: 
    # (csv): 'Mortierella humilis PMI_1414 v1.0'
    # (xml): 'Mortierella humilis PMI_1414 - Glomeribacter phylotype 2 Fungal Standard Draft'
    for element in xml_Assembly.iterdescendants():
        filename = element.attrib["filename"]
        url = element.attrib["url"]
        portal = url.split("/")[2]
        file_size = int(element.attrib["sizeInBytes"])

        if "MitoAssembly" in filename \
            or "MitoScaffolds" in filename \
            or "PrimaryAssemblyScaffolds" in filename \
            or "SecondaryAssemblyScaffolds" in filename \
            or filename in exclude_assemblies \
            or portal in portals_to_remove:
            continue

        try:
            organisms_csv[portal].assembly_file = filename
            organisms_csv[portal].assembly_url = url
            organisms_csv[portal].assembly_size = file_size
        except KeyError:
            print("Warning! Portal {} from xml file not found in csv list".format(portal))
        
        duplicated_names[portal].append(filename)
        
    if len(duplicated_names) != sum([len(x) for x in duplicated_names.values()]):
        print("Warning: the following portals have more than one assembly file!")
        for portal, listfilenames in duplicated_names.items():
            if len(listfilenames) != 1:
                print("{}\t{}".format(portal, ", ".join(listfilenames)))
        print("")

    return


# These files have the "GeneCatalog" magic string but are unwanted
ignore_gff_filenames = {"Aciri1_meta_GeneCatalog_genes_20111216.gff.gz", 
                        "Exoaq1_GeneCatalog_20160901.gff3.gz",
                        "Exoaq1_GeneCatalog_20160828.gff3.gz",
                        "Fonpe1_GeneCatalog_20160901.gff3.gz",
                        "Copmic2_FM1_removed_alleles.gff.gz"}
def get_hardcoded_gffs():
    hgfs = dict()
    tsv_file = Path(__file__).parent / "hardcoded_gff_files.tsv"
    if not tsv_file.is_file():
        print("Missing 'hardcoded_gff_files.tsv'")
        return hgfs
    
    with open(tsv_file) as f:
        for line in f:
            if line[0] == "#":
                continue
            portal, gff_file = line.strip().split("\t")
            hgfs[portal] = gff_file
            
    return hgfs

    
def annotate_gff(outputfolder, organisms_csv, xml_gff, hardcoded_gff_files):
    gff_filenames = defaultdict(list) # each element a pair of filename, datetimeobject
    skipped_filenames = set()

    # If we're to compare the different gff files to get the newest
    # we must convert the timestamp to something we can compare.
    # Because the second parameter in strptime (the format) does not 
    # work with all possible timezones, we must convert to UTC time and 
    # use the %z flag instead of %Z
    UStimezones = {"EST": "-0500", "CST": "-0600", "MST": "-0700", 
        "PST": "-0800", "EDT": "-0400", "CDT": "-0500", "MDT": "-0600", 
        "PDT":"-0700"}
    # Doing this as well for weekdays and months to avoid changing locale setting
    USdaynums = {"Sun":"0", "Mon":"1", "Tue":"2", "Wed":"3", "Thu":"4", 
        "Fri":"5", "Sat":"6"}
    USmonthnums = {"Jan":"01", "Feb":"02", "Mar":"03", "Apr":"04", "May":"05", 
        "Jun":"06", "Jul":"07", "Aug":"08", "Sep":"09", "Oct":"10", 
        "Nov":"11", "Dec":"12"}

    for element in xml_gff.iterdescendants():
        filename = element.attrib["filename"]
        url = element.attrib["url"]
        portal = url.split("/")[2]
        timestamp = element.attrib["timestamp"]
        file_size = int(element.attrib["sizeInBytes"])
        
        if portal in portals_to_remove:
            continue
        
        org = organisms_csv[portal]
        
        # Convert the timestamp-date into a datetime object
        # e.g. 'Sun Oct 12 11:02:03 PDT 2014'
        wday, month, day, HMS, tz, y = timestamp.split(" ")
        stringtimestamp = "{} {} {} {} {} {}".format(USdaynums[wday], 
            USmonthnums[month], day, HMS, UStimezones[tz], y)
        dt_timestamp = datetime.strptime(stringtimestamp, '%w %m %d %H:%M:%S %z %Y')
        
        gff_filenames[portal].append((filename, dt_timestamp))
        
        if portal in hardcoded_gff_files:
            if hardcoded_gff_files[portal] == filename:
                org.gff_file = filename
                org.gff_url = url
                org.gff_timestamp = dt_timestamp
                org.gff_size = file_size
                continue
            else:
                skipped_filenames.add(filename)
                continue

        if filename in ignore_gff_filenames:
            skipped_filenames.add(filename)
            continue
        
        if "proteins" in filename.lower() or "secondary_alleles" in filename.lower():
            # based on example of Aspzo1, these seem to be large files
            # with a compilation of several gff3 files (or secondary alleles)
            skipped_filenames.add(filename)
            continue
        
        # what are these...?
        if filename[-6:] == "gtf.gz" \
                or filename[-3:] == "tgz" \
                or filename[-2:] != "gz":
            skipped_filenames.add(filename)
            continue

        
        if org.gff_file == "":
            # First filename for this portal, record
            org.gff_file = filename
            org.gff_url = url
            org.gff_timestamp = dt_timestamp
            org.gff_size = file_size
        else:
            if filename[-7:] == "gff3.gz" and org.gff_file[-6:] == "gff.gz":
                # had gff, now got gff3. Update
                org.gff_file = filename
                org.gff_url = url
                org.gff_timestamp = dt_timestamp
                org.gff_size = file_size
            elif filename[-7:] == "gff3.gz" and org.gff_file[-7:] == "gff3.gz" or \
                filename[-6:] == "gff.gz" and org.gff_file[-6:] == "gff.gz":
                if dt_timestamp > org.gff_timestamp:
                    # have a file with the same extension. Keep the newest
                    org.gff_file = filename
                    org.gff_url = url
                    org.gff_timestamp = dt_timestamp
                    org.gff_size = file_size
            elif filename[-6:] == "gff.gz" and org.gff_file[-7:] == "gff3.gz":
                # a gff file appeared but we had gff3, skip
                continue
            else:
                print(" annotate_gff unexpected case ({})".format(portal))
                print(org.gff_file)
                print(filename)
                print()
                continue

    repeated_file = o / "List_gene_gff_filenames.txt"
    with codecs.open(repeated_file, "w", encoding='utf-8') as f:
        for portal in sorted(gff_filenames):
            #f.write("{} ({})\n".format(portal, codecs.decode(organisms_csv[portal].name, encoding='utf-8')))
            f.write("{} ({})\n".format(portal, organisms_csv[portal].name))
            f.write("\t{}\n".format(organisms_csv[portal].project_path))
            for filename, dt in gff_filenames[portal]:
                if filename in skipped_filenames:
                    f.write("\t{}\t{} (SKIPPED)".format(dt.strftime("%Y-%m-%d"), filename))
                else:
                    f.write("\t{}\t{}".format(dt.strftime("%Y-%m-%d"), filename))
                    
                if filename == organisms_csv[portal].gff_file:
                    f.write(" *\n")
                else:
                    f.write("\n")
            f.write("\n")
   
    portals_missing_gff_file = [p.portal for p in organisms_csv.values() if p.gff_file == ""]
    if len(portals_missing_gff_file) > 0:
        print("Missing gffs from: {}".format(", ".join(portals_missing_gff_file)))
    

def download_file(url, local_path, cookie):
    """
    Handles downloads for files
    """

    command = ["curl", "https://genome.jgi.doe.gov{}".format(url), "-b", "{}".format(cookie),
               "--retry", "5", "--output", "{}".format(local_path)]
    try:
        sp.run(command, shell=False, stdout=DEVNULL, stderr=STDOUT)
    except CalledProcessError as e:
        print(" Error: {}".format(str(e.output)))
        return False
    else:
        return True


if __name__ == "__main__":
    args = command_parser()

    cookie_path = Path(__file__).parent / "jgi_cookie.txt"

    # Update/download files and exit
    if args.update or args.getgenomelist or args.getxml or args.getprevious:
        if args.update:
            print("Updating NCBI taxonomy database")
            ncbi = NCBITaxa()
            ncbi.update_taxonomy_database()
            print("  Done")
            
        if args.getgenomelist or args.getxml:
            #  Output folder
            o = args.outputfolder
            if not o.is_dir():
                os.makedirs(o, exist_ok=True) # recursive folder creation
        
        if args.getgenomelist:
            print("Downloading genome list (.csv)")
            if get_JGI_genome_list(o):
                print("  Done")
            else:
                print("  Cannot download MycoCosm's genome list...")

        if args.getxml:
            print("Downloading xml file")
            if not JGI_login(cookie_path):
                sys.exit("Error: Cannot log in to JGI...")
            if get_JGI_xml(o, cookie_path):
                print("  Done")
            else:
                print("  Cannot download MycoCosm's xml file...")

        if args.getprevious:
            print("Downloading location list of previously downloaded files")
            p = args.getprevious
            if not p.is_dir():
                sys.exit("Error: given base folder is not a valid folder")
            with open(Path(__file__).parent / "previously_downloaded_files.tsv", "w") as f:
                for gz in p.glob("**/*.gz"):
                    f.write("{}\t{}\n".format(gz.name, gz.parent))
            print("  Done")
        
        sys.exit("All data fetched")
        
    #  Output folder
    o = args.outputfolder
    if not o.is_dir():
        os.makedirs(o, exist_ok=True) # recursive folder creation

    # Read the info from the official list of release genomes
    # It should be the same fungi we'll find defined in the xml
    if not args.csv:
        sys.exit("Error: please download the genome list and use the --csv parameter")
    else:
        csv_path = args.csv
        if not csv_path.is_file():
            sys.exit("Error: cannot open specified csv file")
    
    # The main goal of the script is to populate all the info in each
    # JGI_Project object, then use that info to download/copy the files
    organisms_csv = dict()
    organisms_csv = read_mycocosm_csv(args.csv)
    print("Found {} fungi in the csv genome list\n".format(len(organisms_csv)))

    location_previous = dict()
    if args.previous:
        with open(args.previous) as f:
            for line in f:
                filename, location = line.strip().split("\t")
                location_previous[filename] = location

    if not args.xml:
        sys.exit("Error: please download the xml file and use the --xml parameter")

    if not (args.xml).is_file():
        sys.exit("Error: invalid file specified with --xml")

    # Grab the nodes from the xml that we'll need
    doc = etree.parse(str(args.xml))
    root = doc.getroot()
    for child in root:
        if child.attrib["name"] == "Files":
            Fungi_files = child
            break

    for folder in Fungi_files:
        if folder.attrib["name"] == "Assembly":
            for subfolder in folder:
                if subfolder.attrib["name"] == "Assembled scaffolds (unmasked)":
                    xml_Assembly = subfolder
                    break
    
        if folder.attrib["name"] == "Annotation":
            Annotation_folder = folder
            for subfolder in Annotation_folder:
                if subfolder.attrib["name"] == "Filtered Models (\"best\")":
                    Filtered = subfolder
                    for subsubfolder in Filtered:
                        if subsubfolder.attrib["name"] == "Proteins":
                            xml_Proteins = subsubfolder
                        if subsubfolder.attrib["name"] == "Genes":
                            xml_gff = subsubfolder
    

    # Annotate assembly
    annotate_assembly(organisms_csv, xml_Assembly)
    # Annotate gff files
    hardcoded_gff_files = get_hardcoded_gffs()
    annotate_gff(o, organisms_csv, xml_gff, hardcoded_gff_files)

    # Download (or copy) all data
    if args.simulate:
        print("\nBeginning simulation")
    else:
        print("\nBeginning file download")
        if not JGI_login(cookie_path):
            sys.exit("Error: Cannot log in to JGI...")

    # Make a list of all files for each portal. If there's any error
    # we can start here
    download_list_file = "JGI_download_list_{}.txt".format(time.strftime("%Y-%m-%d", time.localtime()))
    with open(o/download_list_file, "w", encoding='utf-8') as f:
        for portal in organisms_csv:
            fungus = organisms_csv[portal]

            #f.write("{} ({})\n".format(portal, codecs.decode(fungus.name, encoding='utf-8')))
            f.write("{} ({})\n".format(portal, fungus.name))
            f.write("Assembly:\t{}\n".format(fungus.assembly_file))
            f.write("GFF:\t\t{}\n\n".format(fungus.gff_file))
            
    # Finally, get the files
    needed = 0
    copied = 0
    downloaded = 0
    pre_existing = 0
    with open(o / "JGI_taxonomy.tsv", "w", encoding='utf-8') as tf:
        tf.write("Short name\tAccession\tTaxId\tName\tPath\tAssembly file\tGFF file\tlineage\n")

        for portal in organisms_csv:
            fungus = organisms_csv[portal]

            portal = fungus.portal

            if fungus.assembly_file == "" or fungus.gff_file == "":
                print(" Warning! Missing file for {}. See {}".format(portal, download_list_file))
                continue

            needed += 2

            base_folder = Path(fungus.project_path) / portal
            output_folder = o / base_folder
            if not output_folder.is_dir():
                os.makedirs(output_folder, exist_ok=True)

            tf.write("{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format(portal, portal, \
                #fungus.TaxId, codecs.decode(fungus.name, encoding='utf-8'), base_folder, fungus.assembly_file, \
                #fungus.gff_file, fungus.lineage_list))
                fungus.TaxId, fungus.name, base_folder, fungus.assembly_file, \
                fungus.gff_file, fungus.lineage_list))
            if args.simulate:
                continue

            asm = output_folder / fungus.assembly_file
            gff = output_folder / fungus.gff_file

            # Assembly
            if asm.is_file() and asm.stat().st_size > 0.9*fungus.assembly_size:
                # File already there, skipping
                pre_existing += 1
            else:
                if asm in location_previous:
                    try:
                        shutil.copy(location_previous[asm], output_folder)
                    except:
                        sys.exit("Error: cannot copy {} into {}".format(location_previous[asm], output_folder))
                    else:
                        copied += 1
                else:
                    if download_file(fungus.assembly_url, asm, cookie_path):
                        downloaded += 1
                    else:
                        print("Warning: could not download assembly {}".format(asm))

            # GFF
            if gff.is_file() and gff.stat().st_size > 0.9*fungus.gff_size:
                pre_existing += 1
            else: 
                if gff in location_previous:
                    try:
                        shutil.copy(location_previous[gff], output_folder)
                    except:
                        sys.exit("Error: cannot copy {} into {}".format(location_previous[gff], output_folder))
                    else:
                        copied += 1
                else:
                    if download_file(fungus.gff_url, gff, cookie_path):
                        downloaded += 1
                    else:
                        print("Warning: could not download gff {}".format(gff))

            

    print("...done\n")

    print("Portals: {}".format(len(organisms_csv)))
    print("Files needed: {}".format(needed))
    total_got = copied + downloaded + pre_existing
    print(f"Got {copied} (copied) + {downloaded} (downloaded) + {pre_existing} (pre-existing) = {total_got}")

