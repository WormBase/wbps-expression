from Bio import Entrez
from argparse import ArgumentParser

def get_args():
	parser = ArgumentParser()

	parser.add_argument("bioproject", help = "Accession ID for BioProject to extract sample metadata for.")
	parser.add_argument("output_file", help = "Output .json/.tsv file to write sample metadat to.")
	parser.add_argument("-e", "--entrez_email", default = "Matej.Vucak@glasgow.ac.uk", help = "Email address to access Entrez with.")

	args = parser.parse_args()
	return args 

def main():
	args = get_args()

	Entrez.email = args.entrez_email

	search_handle = Entrez.esearch(db="bioproject", retmax = 1, term=args.bioproject, idtype="acc")
	bioproject_record = Entrez.read(search_handle)
	search_handle.close()
	bioproject_id = bioproject_record["IdList"][0]


	elink_handle =  Entrez.elink(dbfrom = "bioproject", db = "biosample", id = bioproject_id)

	record = Entrez.read(elink_handle)
	elink_handle.close()
	biosample_ids = []
	for biosample_link in record[0]["LinkSetDb"][0]["Link"]:
		biosample_id = biosample_link["Id"]
		biosample_ids.append(biosample_id)

	efetch_handle = Entrez.efetch(db="biosample", id=biosample_ids, retmode='xml')
	print(efetch_handle.url)
	biosample_records = Entrez.parse(efetch_handle, validate = False)

	for rec in biosample_records:
		print(rec["Title"])
	efetch_handle.close()



if __name__ == "__main__":
	main()