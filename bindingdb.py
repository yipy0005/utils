import csv
import warnings
from datetime import date
from pathlib import Path
from typing import Optional
from urllib.parse import quote

import requests
import xmltodict
from bs4 import BeautifulSoup
from tap import Tap

warnings.filterwarnings("ignore")

XML_HEADER = '<?xml version="1.0" encoding="UTF-8"?>'

TODAY = date.today().strftime("%y%m%d")


class SimpleArgumentParser(Tap):
    uniprot: Optional[str] = None  # Uniprot ID
    ic50_cutoff: Optional[float] = None  # IC50 (in nM) Cutoff
    smiles: Optional[str] = None  # SMILES string
    similarity_cutoff: Optional[float] = None  # Molecular Similarity Cutoff


def get_ligands_by_uniprot(uniprot: str, ic50_cutoff: float) -> str:
    """
    Queries the BindingDB RESTful API for ligands associated with a given UniProt
    identifier.

    Args:
        uniprot (str): A UniProt identifier.
        ic50_cutoff (float): The IC50 affinity cutoff for the ligands.

    Returns:
        None

    Raises:
        None

    Example Usage:
    get_ligands_by_uniprot("P35355", 100)
    """
    url = (
        "https://bindingdb.org/axis2/services/BDBService"
        f"/getLigandsByUniprot?uniprot={uniprot};{ic50_cutoff}"
    )
    response = requests.get(url)

    xml_content = BeautifulSoup(f"{XML_HEADER}\n{response.text}", "lxml").prettify()
    my_dict = xmltodict.parse(xml_content)

    with open(f"{TODAY}_{uniprot}_{ic50_cutoff}nM.csv", "w") as output_file:
        header = ["SMILES", "Affinity Type", "Affinity (nM)"]
        csv_writer = csv.writer(output_file)
        csv_writer.writerow(header)

        for ligand in my_dict["html"]["body"]["bdb:getligandsbyuniprotresponse"][
            "bdb:affinities"
        ]:
            result = [
                ligand["bdb:smiles"],
                ligand["bdb:affinity_type"],
                ligand["bdb:affinity"],
            ]

            csv_writer.writerow(result)

    return (
        "Results have been saved to "
        f"{Path(f'{TODAY}_{uniprot}_{ic50_cutoff}nM.csv').parent.absolute()}"
        f"/{TODAY}_{uniprot}_{ic50_cutoff}nM.csv!"
    )


def get_target_by_compound(smiles: str, similarity_cutoff: float):
    """
    Queries the BindingDB RESTful API for compounds similar to a given SMILES string and
    their protein targets.

    Args:
        smiles (str): The SMILES string for the query compound.
        similarity_cutoff (float): The similarity cutoff for finding similar compounds.

    Returns:
        None

    Raises:
        None

    Example Usage:
    get_target_by_compound("CCC[Ni+](C)(C)CCn1nncc1COc1cc(=O)n(C)c2ccccc12", 0.85)
    """
    url = (
        "https://bindingdb.org/axis2/services/BDBService"
        f"/getTargetByCompound?smiles={quote(smiles)}&cutoff={similarity_cutoff}"
    )
    response = requests.get(url)

    xml_content = BeautifulSoup(f"{XML_HEADER}\n{response.text}", "lxml").prettify()
    my_dict = xmltodict.parse(xml_content)

    with open(
        f"{TODAY}_targetSearchByCpd_{similarity_cutoff}nM.csv", "w"
    ) as output_file:
        header = [
            "Target [Species]",
            "SMILES (Tanimoto Similarity)",
            "Affinity Type",
            "Affinity (nM)",
        ]
        csv_writer = csv.writer(output_file)
        csv_writer.writerow(header)
        for target in my_dict["html"]["body"]["bdb:gettargetbycompoundresponse"][
            "bdb:affinities"
        ]:
            result = [
                f"{target['bdb:target']} [{target['bdb:species']}]",
                f"{target['bdb:smiles']} ({target['bdb:tanimoto']})",
                target["bdb:affinity_type"],
                target["bdb:affinity"],
            ]

            csv_writer.writerow(result)

        return (
            "Results have been saved to "
            f"{Path(f'{TODAY}_{smiles}_{similarity_cutoff}nM.csv').parent.absolute()}"
            f"/{TODAY}_targetSearchByCpd_{similarity_cutoff}nM.csv!"
        )


def main(
    uniprot: Optional[str],
    ic50_cutoff: Optional[float],
    smiles: Optional[str],
    similarity_cutoff: Optional[float],
):
    if uniprot and ic50_cutoff:
        print(get_ligands_by_uniprot(uniprot, ic50_cutoff))
    elif smiles and similarity_cutoff:
        print(get_target_by_compound(smiles, similarity_cutoff))


if __name__ == "__main__":
    args = SimpleArgumentParser().parse_args()
    main(
        args.uniprot,
        args.ic50_cutoff,
        args.smiles,
        args.similarity_cutoff,
    )
