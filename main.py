from blast import *
from input_parser import parse_cli
from input_parser import ParseInput
from output import write_output
from output import create_images
def main():
    parse_res = parse_cli()
    input_file = parse_res["input_file_name"]
    output_file = parse_res["output_file_name"]
    pichia_genes = parse_res["list_genes"]

    print("============ Start ============")

    if len(pichia_genes) == 0:
        pichia_genes = ParseInput(input_file=input_file).run()
    else:
        if any((not gene.startswith("PAS") for gene in pichia_genes)):
            raise ValueError("Гены должны начинаться на PAS")

    result = Blast(pichia_genes=pichia_genes).protein_blast()
    write_output(result)
    create_images()

    print("============= Done ============")

if __name__ == "__main__":
    main()