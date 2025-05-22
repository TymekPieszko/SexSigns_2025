import pytest
import newick
from sexsigns_functions.calc import rename_nodes, get_tree_category

# Create inputs
# Execute code
# Compare output with expected result


class Unbalanced4TipNewick:
    def __init__(self, A, B, C, D):
        self.A = A
        self.B = B
        self.C = C
        self.D = D

    def get_string(self):
        return f"((({self.A}:0.1,{self.B}:0.1):0.1,{self.C}:0.1):0.1,{self.D}:0.1);"


class Balanced4TipNewick:
    def __init__(self, A, B, C, D):
        self.A = A
        self.B = B
        self.C = C
        self.D = D

    def get_string(self):
        return f"(({self.A}:0.1,{self.B}:0.1):0.1,({self.C}:0.1,{self.D}:0.1):0.1);"


@pytest.fixture
def unbalanced_tree_before_renaming():
    return Unbalanced4TipNewick("n123", "n156", "n124", "n157")


@pytest.fixture
def balanced_tree_before_renaming():
    return Balanced4TipNewick("n123", "n156", "n124", "n157")


@pytest.fixture
def u0123():
    return Unbalanced4TipNewick("0", "1", "2", "3")


class TestRenameNodes:
    def test_unbalanced_tree(self, unbalanced_tree_before_renaming):
        target_names = ["0", "1", "2", "3"]
        tree = newick.loads(unbalanced_tree_before_renaming.get_string())[0]
        tree = rename_nodes(tree)
        new_names = [node.name for node in tree.walk() if node.name is not None]
        new_names.sort()
        assert (
            new_names == target_names
        ), f"Nodes should be renamed to {target_names}, got {new_names}"

    def test_balanced_tree(self, balanced_tree_before_renaming):
        target_names = ["0", "1", "2", "3"]
        tree = newick.loads(balanced_tree_before_renaming.get_string())[0]
        tree = rename_nodes(tree)
        new_names = [node.name for node in tree.walk() if node.name is not None]
        new_names.sort()
        assert (
            new_names == target_names
        ), f"Nodes should be renamed to {target_names}, got {new_names}"

# class Balanced4TipNewick:
#     def __init__(self, A, B, C, D, br_means, br_sds):
#         self.A = A
#         self.B = B
#         self.C = C
#         self.D = D
#         self.br_means = br_means
#         self.br_sds = br_sds

#     def get_string(self):
#         br_means = self.br_means
#         np.random.shuffle(br_means)
#         # print(mean_lengths)
#         m1 = br_means[0]
#         m2 = br_means[1]
#         s1 = self.br_sds[0]
#         s2 = self.br_sds[1]
#         return f"(({self.A}:{norm(m1, s1)},{self.B}:{norm(m1, s1)}):0.1,({self.C}:{norm(m2, s2)},{self.D}:{norm(m2, s2)}):0.1);"