from polars import col, DataFrame

class Fastg_writer():

    @classmethod
    def parse_main_descriptor(cls, row):

        name      = row[0][:len(row[0])-1]
        length    = str(row[1])
        coverage  = f"{row[2]:.6f}"
        rc        = "" if "+" in row[0] else "'"

        return f"EDGE_{name}_length_{length}_cov_{coverage}{rc}"
    
    @classmethod
    def parse_descriptor(cls, edge_data, row):
        """
        Creates an entry description in .fastg-format
        for a row in the edge_data dataframe.

        Parameters
        ----------
        edge_data: polars.DataFrame
            A dataframe of edge data.
            A[key][i] corresponds to attribute "key" of the i-th edge of an assembly graph
        file_path: tuple
            A row in edge_data
        mode: str
            
        """
        descriptor =  ">" + cls.parse_main_descriptor(row)
        neighbors  = row[4]

        if len(neighbors) > 0:

            neighbor_rows        = edge_data.filter(col('name').is_in(neighbors))
            neighbor_descriptors = [cls.parse_main_descriptor(neighbor_row) for neighbor_row in neighbor_rows.iter_rows()]
            descriptor          += ":" + ",".join(neighbor_descriptors)

        return descriptor + ";"
        
    @classmethod
    def write(cls, edge_data: DataFrame, file_path: str, mode='a')->None:
        """
        Writes assembly graph from dataframe to a file

        Parameters
        ----------
        edge_data: polars.DataFrame
            A dataframe of edge data.
            A[key][i] corresponds to attribute "key" of
            the i-th edge of an assembly graph.
        file_path: str
            File to write to
        mode: str (optional)
            'a' for appending (default)
            'w' for overwriting
        """
        with open(file_path, mode) as file:
            for row in edge_data.iter_rows():

                descriptor = cls.parse_descriptor(edge_data, row)
                sequence   = row[3]

                file.write(sequence+"\n")
                file.write(descriptor+"\n")