import pickle
import hashlib
import os

nease_folder = '/home/elias/DIGGER/container/nease_events/'
files = os.listdir(nease_folder)
example_data = pickle.load(open(os.path.join(nease_folder + files[0]), 'rb'))


byte_stream = pickle.dumps(example_data)
hash = hashlib.sha256(byte_stream).hexdigest()


def compute_file_hash(file_path):
    hash_object = hashlib.sha256()
    with open(file_path, 'rb') as file:
        # Read the file in chunks to avoid memory issues with large files
        for chunk in iter(lambda: file.read(4096), b""):
            hash_object.update(chunk)
    return hash_object.hexdigest()


for filename in os.listdir(nease_folder):
    file_path = os.path.join(nease_folder, filename)
    if os.path.isfile(file_path):
        existing_file_hash = compute_file_hash(file_path)
        if existing_file_hash == hash:
            print(f"Found identical file: {filename}")
            break
        else:
            # # If no identical file is found, write the new file
            # new_file_path = os.path.join(nease_folder, 'new_file.pkl')
            # with open(new_file_path, 'wb') as new_file:
            # new_file.write(byte_stream)
            print(f"New file write")