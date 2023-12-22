# test_sequence_properties.py

import pytest
from ReadGenbank import ReadGenbank
import json

def test_circular_sequence_p20N31():
    reader = ReadGenbank()
    result_json = reader.run("p20N31.ape", [])
    result = json.loads(result_json)
    assert result["circular_linear"] == "circular"
    assert result["length"] == 2911

def test_linear_sequence_SCU49845():
    reader = ReadGenbank()
    result_json = reader.run("SCU49845.gb", [])
    result = json.loads(result_json)
    assert result["circular_linear"] == "linear"
    assert result["length"] == 5028

# Additional tests can be added for other sequences and properties
#  Add one from benchling