# Some utility classes to represent a PDB structure


class Atom:
    """
    A simple class for an amino acid residue
    """

    def __init__(self, type_var):
        self.type = type_var
        self.coords = (0.0, 0.0, 0.0)

    # Overload the __repr__ operator to make printing simpler.
    def __repr__(self):
        return self.type


class Residue:
    """
    A simple class for an amino acid residue
    """

    def __init__(self, type_var, number):
        self.type = type_var
        self.number = number
        self.atoms = []

    # Overload the __repr__ operator to make printing simpler.
    def __repr__(self):
        return self.type + " " + self.number


class ActiveSite:
    """
    A simple class for an active site
    """

    def __init__(self, name):
        self.name = name
        self.residues = []

    # Overload the __repr__ operator to make printing simpler.
    def __repr__(self):
        return self.name
