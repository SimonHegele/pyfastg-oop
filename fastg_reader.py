from re     import compile
from polars import DataFrame

class FileError(Exception):
    pass

class Fastg_reader():
    """
    Reads .fastg-files and converts them into networkx.DiGraph()-objects
    """

    edge_descriptor_pattern = compile(r"^EDGE_(?P<name>[a-zA-Z\d]+)_length_(?P<length>\d+)_cov_(?P<coverage>[\d\.]+)$")

    @classmethod
    def check_file_ending(cls, fastg_file: str):

        if not (fastg_file.endswith('.fastg')):
            raise FileError('The path provided does not end with .fastg.\nIf it is a .fastg file please rename it and try again')
    
    @classmethod
    def read_file(cls, fastg_file: str)->tuple:
        '''
        Parameters:
        -----------
        fastg_file: str
            path to .fastg-file
        
        Returns:
        --------
        list:
            A list that contains two-tuples of lines from the .fastg-file,
            each corresponding to an entry of the file
        list:
            A list that containes indices referring to lines in the .fastg-file

        '''
        with open(fastg_file, "r") as file:
            # reading file content
            lines   = file.readlines()
            lines   = [line.strip() for line in lines]
            # splitting file content into individual entries
            indices = [i for i, line in enumerate(lines) if line.startswith('>')] + [len(lines)]
            entries = [(lines[indices[i]], "".join([lines[j] for j in range(indices[i]+1,indices[i+1])])) for i in range(len(indices)-1)]
            return entries, indices
        
    @classmethod
    def check_entry_format(cls, entry: tuple, i: int)->None:
        """
        Checks if the format 
        Parameters:
        -----------
        entry: tuple
            a two-tuple of strings representing an edge of the assembly graph where
            entry[0] is an edge descriptor in SPAdes .fastg-format
            entry[1] the nucleotide sequence of the edge
        i: int
            Line in the file corresponding to the entry

        Raises:
        -------
        FileError
            Reports various potential formatting issues of entries within a .fastg-file
        """
        msg = None
        if "[" in entry[0]:
            msg = f"Line {i}: Descriptor: [] notation not supported" + "\n" + entry[0]
        if not entry[0].endswith(";"):
            msg = f"Line {i}: Descriptor: must end with ;" + "\n" + entry[0]
        if not entry[0].endswith(";"):
            msg = f"Line {i}: Descriptor: not more than one : allowed" + "\n" + entry[0]
        if len([char for char in entry[1] if not char in "ACGTU"])>0:
            msg = f"Line {i+1} Sequence: not a valid nucleotide sequence (contains prohibited characters)" + "\n" + entry[0] + "\n" + entry[1]
        if msg is not None:
            raise FileError(msg)
    
    @classmethod
    def extract_edge_properties(cls, edge_descriptor)->dict:
        """
        Returns attributes of a .fastg edge descriptor as a dictionary.

        Example:
            In:  extract_node_attrs("EDGE_3_length_100_cov_28.087'")
            Out: {"name": "3-", "length": 100, "coverage": 28.087}

        Parameters
        ----------
        declaration: str
            edge descriptor

        Returns
        -------
        dict
            {"name": str, "length": int, "coverage": float}
            edge properties
        """
        rc = True if edge_descriptor.endswith("'") else False # Probably refers to the strand of the strand
        if rc:
            edge_descriptor = edge_descriptor[:-1]

        match   = cls.edge_descriptor_pattern.search(edge_descriptor)
        if match is None:
            msg = f'File invalid. Check your descriptors, only SPAdes dialect is supported {edge_descriptor}'
            raise FileError(msg)
        
        name     = match.group("name")
        name     = name+"-" if rc else name+"+"
        length   = int(match.group("length"))
        coverage = float(match.group("coverage"))

        return {"name": name, "length": length, "coverage": coverage}

    @classmethod
    def entry_to_dict(cls, entry: tuple)->dict:
        """
        For an edge of an assembly graph in .fastq format a corresponding data is extracted
        and returned as a dictionary

        Remark: I will omit the consistency check that was included in the original pyfastg 
                and not because I'm too lazy to rewrite it. Let's trust the SPAdes authors,
                for now.

        Parameters:
        -----------
        entry: tuple
            a two-tuple of strings representing an edge of the assembly graph where
            entry[0] is an edge descriptor in SPAdes .fastg-format
            entry[1] the nucleotide sequence of the edge

        Returns:
        --------
        dict:
            Representation of an edge and its attributes as dictionary
            {"name": str, "length": int, "coverage": float, neighbors: list}
        """
        edge_descriptor = entry[0].rstrip(";")
        edge_properties = cls.extract_edge_properties(edge_descriptor.split(":")[0][1:])
        if ":" not in edge_descriptor:
            # orphaned or terminal edge
            neighbors = []
        else:
            # has at least one adjacent edge
            neighbors = edge_descriptor.split(":")[1].split(",")
            try:
                neighbors = [cls.extract_edge_properties(edge_descriptor)["name"] for edge_descriptor in neighbors]
            except:
                print(entry[0])
                print(entry[1])
                x = 0/0
        edge_properties.update({"sequence": entry[1], "neighbors": neighbors})

        return edge_properties
   
    def parse_fastg(self, fastg_file: str):
        '''
        Parameters:
        -----------
        fastg_file: str
            path to .fastg-file
        '''
        self.check_file_ending(fastg_file)
        
        entries, indices = self.read_file(fastg_file)

        node_data   = {
            "name":       [None for _ in entries],
            "length":     [None for _ in entries],
            "coverage":   [None for _ in entries],
            "sequence":   [None for _ in entries],
            "neighbors":  [None for _ in entries],
            }

        for i, entry in enumerate(entries):

            self.check_entry_format(entry, indices[i])
            entry_dict = self.entry_to_dict(entry)
            for key in node_data.keys():
                node_data[key][i] = entry_dict[key]
            
        return DataFrame(node_data)