import os

# Determine project root = parent directory of scripts/
project_root = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))

def get_path(*parts):
    """
    Build an absolute path relative to the project root.

    Example:
        get_path("data", "BAC88139.1_top150.fasta")
    returns:
        /home/.../cyano_rhodopsins/data/BAC88139.1_top150.fasta
    """
    return os.path.join(project_root, *parts)
