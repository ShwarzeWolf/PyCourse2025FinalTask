import pandas as pd
import rdkit
from rdkit import Chem


class SmilesFileParsingException(Exception):
    pass


def check_file_corections(file_name):
    df = pd.read_csv(file_name)
    columns = {column.lower().strip() for column in df.columns}
    # rename columns
    if "smiles" not in columns:
        raise SmilesFileParsingException("No SMILES column in the file")

    if "molecule name" not in columns:
        raise SmilesFileParsingException("No Molecule Name column in the file")


def cheak_content_corections(file_name):
    df = pd.read_csv(file_name)
    df["molecule"] = df["SMILES"].apply(
        lambda smiles: Chem.MolFromSmiles(smiles)
    )
    valid_molecules_df = df[df["molecule"].notna()]
    not_valid_molecules_df = df[df["molecule"].isna()]

    valid_molecules_df = valid_molecules_df.drop(["molecule"], axis=1)
    valid_molecules_df.drop_duplicates(inplace=True)
    valid_molecules_df.to_csv("data_03_2024_filtered.csv", index=False)

    return not_valid_molecules_df
    
