from abc import ABC, abstractmethod

class FastaReader(ABC):
    """
    Reading .fasta-files
    """
    
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
            A list that contains two-tuples of lines from the .fasta-file,
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

        @abstractmethod
        def parse(self, file):
            pass
