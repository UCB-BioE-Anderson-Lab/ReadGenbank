import pytest
from ReadGenbank import ReadGenbank
import json

def test_full_sequence():
    reader = ReadGenbank()
    result_json = reader.run("SCU49845.gb", ["full_sequence"], [0, 30])
    result = json.loads(result_json)  # Parse JSON string into a Python dictionary

    expected_sequence = "GATCCTCCATATACAACGGTATCTCCACCT"  # Replace with the actual expected sequence

    # Assert that the 'full_sequence' key is in the result and equals the expected sequence
    assert "full_sequence" in result
    assert result["full_sequence"] == expected_sequence

# def test_features():
#     reader = ReadGenbank()
#     result_json = reader.run("SCU49845.gb", ["features"], [600, 700])
#     result = json.loads(result_json)  # Parse JSON string into a Python dictionary

#     # Check if features are in the result
#     assert "features" in result
#     features = result["features"]

#     # Check if the specific feature is present in the extracted features
#     found_feature = False
#     for feature in features:
#         if feature['type'] == 'CDS' and feature['location'] == "87-700" and feature['qualifiers'].get('gene') == ['AXL2']:
#             found_feature = True
#             break

#     assert found_feature, "CDS feature for AXL2 gene in range 600-700 not found"

# Additional tests for other features and data types can be added similarly
