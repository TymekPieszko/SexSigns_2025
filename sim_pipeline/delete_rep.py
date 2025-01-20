import os


def delete_rep_files(directory):
    """
    Deletes all files containing 'rep' in their filenames within the specified directory and its subdirectories.

    Args:
        directory (str): The path to the directory where the search for files will begin.

    Raises:
        Exception: If there is an error deleting a file, it will be caught and printed.

    Example:
        To run this script, use the following command:

        ```python
        delete_rep_files('/path/to/your/directory')
        ```
    """
    for root, dirs, files in os.walk(directory):
        for file in files:
            if "rep" in file:
                file_path = os.path.join(root, file)
                try:
                    os.remove(file_path)
                    print(f"Deleted: {file_path}")
                except Exception as e:
                    print(f"Error deleting {file_path}: {e}")


if __name__ == "__main__":
    directory = "/data/zool-barralab/scro4331/chapter1/phasing_bias/simulation_pipeline_NEW/sim_output"  # Replace with your directory path
    delete_rep_files(directory)
