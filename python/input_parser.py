from typing import List, Optional
import pandas as pd
import argparse
import sys

class ParseInput:
    def __init__(self, input_file: Optional[str] = None):
        self.input_file = input_file

    def pandas_parser(self, separator):
        data = pd.read_csv(self.input_file, sep=separator)

        if "id" not in data.columns:
            if "gene" not in data.columns:
                raise ValueError(
                    'В таблице должен быть столбец с названием либо "gene", либо "id"'
                )
            else:
                pichia_genes = data["gene"]
        else:
            pichia_genes = data["id"]

        if pichia_genes.apply(lambda x: not x.startswith("PAS")).any():
            raise ValueError("Названия генов должны начинаться на PAS")

        return pichia_genes.to_list()

    def csv_input(self) -> List[str]:
        return self.pandas_parser(",")

    def tsv_input(self) -> List[str]:
        return self.pandas_parser("\t")

    def run(self) -> List[str]:
        if self.input_file is None:
            raise ValueError("Нет файла с данными")

        if self.input_file.endswith("csv"):
            return self.csv_input()

        if self.input_file.endswith("tsv"):
            return self.tsv_input()

        raise ValueError("Формат файла должен быть csv или tsv. Либо передайте на вход список генов.")


def parse_cli():
    """
    :return: Возвращает словарь с ключами input_file_name, output_file_name и list_genes.
            input_file_name - входной файл
            output_file_name - файл на запись
            list_genes - список генов
    """
    input_file_name = None
    output_file_name = None
    list_genes = []

    parser = argparse.ArgumentParser(
        description="BLAST Picha pastoris genes (PAS_...) vs Saccharomyces cerevisiae"
    )

    parser.add_argument("args", nargs=argparse.REMAINDER)

    parser.add_argument(
        "stdin", nargs="?", type=argparse.FileType("r"), default=sys.stdin
    )

    parser.add_argument(
        "-f",
        "--file",
        type=argparse.FileType("r"),
        default=sys.stdin,
        help='File with genes. Should be column "id"',
    )

    parser.add_argument(
        "-o",
        "--output",
        type=argparse.FileType("w"),
        default=sys.stdout,
        help="Output with blast result",
    )

    args = parser.parse_args().args
    stdin = parser.parse_args().stdin
    input_file = parser.parse_args().file
    output_file = parser.parse_args().output

    if input_file.name == "<stdin>":
        # Нет файла на вход
        if not sys.stdin.isatty():
            stdin = stdin.read().splitlines()
        else:
            stdin = []

        list_genes = args + stdin
    else:
        input_file_name = input_file.name

    if output_file.name != "<stdout>":
        output_file_name = output_file.name

    return {
        "input_file_name": input_file_name,
        "output_file_name": output_file_name,
        "list_genes": list_genes,
    }
