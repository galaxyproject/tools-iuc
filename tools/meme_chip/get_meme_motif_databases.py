import os


def get_meme_motif_database_options(file_path):
    options = []
    if not os.path.isdir(file_path):
        return options
    for i, file_name in enumerate(os.listdir(file_path)):
        full_path = os.path.join(file_path, file_name)
        options.append((file_name, full_path, i == 0))
    return options
