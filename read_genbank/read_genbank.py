from Bio import SeqIO
import json

class ReadGenbank:
    """
    Parses a GenBank file and returns specified data as a JSON string.
    """
    def run(self, file_path, data_types=None, sequence_range=None):
        """
        Parses a GenBank file and returns specified data based on data types and sequence range.

        Parameters:
        file_path (str): Path to the GenBank file. Supported extensions: .str, .ape, .gb, .seq.
        data_types (list of str): Data types to extract (e.g., 'full_sequence', 'features').
        sequence_range (tuple of int, optional): Start and end positions of the sequence range.

        Returns:
        str: JSON string containing requested data.
            - 'circular_linear': 'circular' or 'linear', indicating the topology of the sequence.
            - 'length': Total length of the sequence.
            - 'full_sequence': The specified range of the DNA sequence.
            - 'features': List of features within the specified range, each as a dictionary.

        Raises:
        ValueError: If file extension is not supported or range exceeds sequence length.
        """
        if data_types is None:
            data_types = []  # Set a default value if none provided


        # Validate file extension
        if not file_path.endswith(('.str', '.ape', '.gb', '.seq')):
            raise ValueError("Invalid file extension. Supported extensions are .str, .ape, .gb, and .seq.")

        # Read the file and parse it
        try:
            with open(file_path, "r") as file:
                record = SeqIO.read(file, "genbank")
        except Exception as e:
            return json.dumps({"error": f"Error reading file: {e}"})

        # Validate the sequence range
        if sequence_range:
            start, end = sequence_range
            if start < 0 or end > len(record.seq):
                raise ValueError("Range exceeds the length of sequence")

        results = {"circular_linear": "circular" if record.annotations.get("topology", "linear") == "circular" else "linear",
                   "length": len(record.seq)}

        # Full sequence retrieval
        if "full_sequence" in data_types:
            if sequence_range is None:
                start, end = 0, len(record.seq)
            sequence = str(record.seq[start:end])
            results["full_sequence"] = sequence

        # Features retrieval
        if "features" in data_types:
            features = []
            for feature in record.features:
                if sequence_range:
                    feature_start = int(feature.location.start)
                    feature_end = int(feature.location.end)
                    if (feature_start < sequence_range[1] and feature_end > sequence_range[0]):
                        partial = feature_start < sequence_range[0] or feature_end > sequence_range[1]
                        relative_start = max(feature_start, sequence_range[0]) - sequence_range[0]
                        relative_end = min(feature_end, sequence_range[1]) - sequence_range[0]
                        feature_info = {
                            'type': feature.type,
                            'location': f"{relative_start}-{relative_end}",
                            'partial': partial,
                            'qualifiers': feature.qualifiers
                        }
                        features.append(feature_info)
                else:
                    features.append({
                        'type': feature.type,
                        'location': str(feature.location),
                        'qualifiers': feature.qualifiers
                    })

            results["features"] = features

        return json.dumps(results, indent=4)

# Example usage
reader = ReadGenbank()
result = reader.run("seqs/p20N31.ape", ["full_sequence"])
print(result)
