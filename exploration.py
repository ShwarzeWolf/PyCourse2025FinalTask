# pip install pandas

import pandas as pd
import rdkit
from rdkit import Chem


class SmilesFileParsingException(Exception):
    pass


# numpy pandas
# параллельное выполнение операций
df = pd.read_csv("data_03_2024.csv")
columns = {column.lower().strip() for column in df.columns}
# rename columns
if "smiles" not in columns:
    raise SmilesFileParsingException("No SMILES column in the file")

if "molecule name" not in columns:
    raise SmilesFileParsingException("No Molecule Name column in the file")

df["molecule"] = df["SMILES"].apply(
    lambda smiles: Chem.MolFromSmiles(smiles)
)

valid_molecules_df = df[df["molecule"].notna()]
not_valid_molecules_df = df[df["molecule"].isna()]

valid_molecules_df = valid_molecules_df.drop(["molecule"], axis=1)
valid_molecules_df.drop_duplicates(inplace=True)
valid_molecules_df.to_csv("data_03_2024_filtered.csv", index=False)

print(
    f"There are some molecules which were not parsed: "
    f"{not_valid_molecules_df['molecule name'].to_list()}."
    "Those molecules were removed from the file."
    "Please change the value if you want to have it in the search."
)