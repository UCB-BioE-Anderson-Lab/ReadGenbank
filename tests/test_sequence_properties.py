# test_sequence_properties.py

import pytest
from ..read_genbank.ReadGenbank import ReadGenbank
import json

def test_linear_sequence_from_NCBI():
    reader = ReadGenbank()
    result_json = reader.run("seqs/SCU49845.gb", [])
    result = json.loads(result_json)
    assert result["circular_linear"] == "linear"
    assert result["length"] == 5028

def test_circular_sequence_from_ApE():
    reader = ReadGenbank()
    result_json = reader.run("seqs/p20N31.ape", [])
    result = json.loads(result_json)
    assert result["circular_linear"] == "circular"
    assert result["length"] == 2911

def test_linear_sequence_from_benchling():
    reader = ReadGenbank()
    result_json = reader.run("seqs/pbr322_egfr.gb", [])
    result = json.loads(result_json)
    assert result["circular_linear"] == "circular"
    assert result["length"] == 5968

def test_linear_sequence_no_datatype():
    reader = ReadGenbank()
    result_json = reader.run("seqs/SCU49845.gb")
    result = json.loads(result_json)
    assert result["circular_linear"] == "linear"
    assert result["length"] == 5028