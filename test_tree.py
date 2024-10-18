import os


def build_directory_tree(root_dir):
    tree = {}
    for entry in os.listdir(root_dir):
        path = os.path.join(root_dir, entry)
        if os.path.isdir(path):
            tree[entry] = build_directory_tree(path)  # Recursively build tree for subdirectories
        else:
            tree[entry] = None  # Files are leaf nodes
    return tree

# Example Usage:
root_directory = "2024"
directory_tree = build_directory_tree(root_directory)

for x in directory_tree:
    print(x)

print(directory_tree["0912_20K_H2O_32min"]==None)
print(type(directory_tree))