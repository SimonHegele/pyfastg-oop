from re     import compile, finditer
from polars import DataFrame

from reader_fasta import FastaReader

class FileError(Exception):
    pass

class BcalmReader(FastaReader):
    """
    Reading .fasta-files from bcalm and returning it's content as a polars.DataFrame
    """
         
    @classmethod
    def entry_to_dict(cls, entry: tuple)->dict:
        """
        For an entry of a .fasta-file from bcalm the corresponding edge data is extracted
        using regex and returned as a dictionary


        Parameters:
        -----------
        entry: tuple
            a two-tuple of strings representing an edge of the assembly graph where
            entry[0] header (e.g. ">2 LN:i:33 KC:i:231 km:f:77.0   L:-:10:-  L:+:0:- L:+:11:-")
            entry[1] nucleotide sequence

        Returns:
        --------
        dict:
            Representation of an edge and its attributes as dictionary
        """
        
        # dictionary of regex expressions
        r = {
            "name":            compile(r"^>(?P<ID>\d+)"),
            "length":          compile(r"LN:i:(?P<LENGTH>\d+)"),
            "total abundance": compile(r"KC:i:(?P<TOTAL_ABUNDANCE>\d+)"),
            "avg. abundance":  compile(r"km:f:(?P<AVG_ABUNDANCE>\d+\.\d)"),
            "edges":           compile(r"(?P<EDGE>L:[-+]:\d:[-+])"),
            "neighbors":       compile(r"L:[-+]:(?P<NEIGHBOR>\d):[-+]")
            }
        
        return {
            "name":            r["name"].search(entry[0]).group("ID"),
            "length":          r["length"].search(entry[0]).group("LENGTH"),
            "total abundance": int(r["total abundance"].search(entry[0]).group("TOTAL_ABUNDANCE")),
            "avg. abundance":  float(r["avg. abundance"].search(entry[0]).group("AVG_ABUNDANCE")),
            "coverage":        float(r["avg. abundance"].search(entry[0]).group("AVG_ABUNDANCE")),
            "sequence":        entry[1],
            "edges":           [match.group("EDGE") for match in finditer(r["edges"], entry[0])],
            "neighbors":       [match.group("NEIGHBOR") for match in finditer(r["neighbors"], entry[0])]
                }
   
    def parse(self, file: str):
        '''
        Implements abstract method parse inherited from FastaReader

        Parameters:
        -----------
        fastg_file: str
            path to .fastg-file
        '''
        
        entries, indices = self.read_file(file) # Method inherited from FastaReader

        node_data   = {
            "name":            [None for _ in entries],
            "length":          [None for _ in entries],
            "total abundance": [None for _ in entries],
            "avg. abundance":  [None for _ in entries],
            "coverage":        [None for _ in entries],
            "sequence":        [None for _ in entries],
            "edges":           [None for _ in entries],
            "neighbors":       [None for _ in entries],
            }

        for i, entry in enumerate(entries):

            entry_dict = self.entry_to_dict(entry)
            for key in node_data.keys():
                node_data[key][i] = entry_dict[key]
            
        return DataFrame(node_data)