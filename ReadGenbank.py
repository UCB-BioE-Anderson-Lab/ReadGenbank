from Bio import SeqIO
import json

class ReadGenbank:
    def run(self, file_path, data_types, sequence_range=None):
        """
        Parses a GenBank (or similar format) file and returns specified data based on the data types and sequence range.
        """
        # Validate file extension
        if not file_path.endswith(('.str', '.ape', '.gb', '.seq')):
            raise ValueError("Invalid file extension. Supported extensions are .str, .ape, .gb, and .seq.")

        # Read the file and parse it
        try:
            with open(file_path, "r") as file:
                record = SeqIO.read(file, "genbank")
        except Exception as e:
            return f"Error reading file: {e}"

        results = []

        # Full sequence retrieval
        if "full_sequence" in data_types:
            if sequence_range:
                start, end = sequence_range
            else:
                start, end = 0, len(record.seq)
            sequence = str(record.seq[start:end])
            results.append(f"Full Sequence [{start}:{end}]:\n{sequence}\n")

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

            pretty_features = json.dumps(features, indent=4)
            results.append(f"Features:\n{pretty_features}\n")

        # Other data types can be implemented similarly

        return "\n".join(results)

# Example usage
reader = ReadGenbank()

# Example 1: Retrieve full sequence and features within a range
try:
    result = reader.run("/Users/jcaucb/Downloads/trashy/SCU49845.gb", ["full_sequence", "features"], [0, 200])
    print(result)
except ValueError as e:
    print(e)
