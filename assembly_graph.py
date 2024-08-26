from networkx          import DiGraph, weakly_connected_components
from polars            import col, DataFrame
from typing_extensions import Iterator, Self

class Assembly_graph(DiGraph):
    '''
    Networkx representation of assembly graphs.

    Here, the actual edges of the assembly graph are represented as nodes and vice versa
    '''

    def __init__(self, edge_data: DataFrame):
        """
        Initializes the assembly graph

        Parameters
        ----------
        edge_data: polars.DataFrame
            A dataframe of edge data.
            A[key][i] corresponds to attribute "key" of the i-th edge of an assembly graph
        """

        super().__init__()

        self.add_nodes_from(edge_data['name'])

        for i, node in enumerate(edge_data['name']):
            for j, neighbor in enumerate(edge_data["neighbors"][i]):
                self.add_edge(node, neighbor)

        self.assembly_data = edge_data

    def get_assembly_data(self, nodes):
        """
        Parameter:
        nodes: collection
            A subset of the graphs nodes.
            If None all nodes of the graph are chosen

        Returns:
        --------
        dict
            A dictionary of edge data.
            A[key][i] corresponds to property "key" of the i-th edge of an assembly (sub-) graph
        """
        return self.assembly_data.filter(col('name').is_in(nodes))

    def component_graph(self, component: set)->Self:
        """
        Parameters:
        -----------
        component: set
            A set of node names of the assembly graph

        Returns:
        --------
        Assembly_graph
            Subgraph of the assembly graph 
            with all nodes provided in components and corresponding edges.
        """
        return Assembly_graph(self.get_assembly_data(component))

    def component_graphs(self)->Iterator:
        """
        Returns:
        iterator:
            A list of Assembly_graph objects, each corresponding to one of the
            disconnected components of the graph.
        """
        components       = sorted(weakly_connected_components(self), key=len, reverse=True)
        for c in components:
            yield self.component_graph(c)