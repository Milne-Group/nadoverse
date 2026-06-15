import json

import pytest

from nadoverse.registry import REGISTRY, NadoTool, all_tools, get_tool, to_json


def test_registry_loads():
    assert len(REGISTRY) > 0


def test_all_tools_returns_copy():
    tools = all_tools()
    assert tools == REGISTRY
    assert tools is not REGISTRY


def test_all_tools_have_required_fields():
    required = {"name", "description", "pypi_name", "install_extra", "repo_url", "min_python"}
    for tool in REGISTRY:
        missing = required - set(vars(tool))
        assert not missing, f"{tool.name} missing fields: {missing}"


def test_cli_tools_have_command():
    cli_tools = [t for t in REGISTRY if t.cli_command is not None]
    assert len(cli_tools) > 0
    for tool in cli_tools:
        assert isinstance(tool.cli_command, str)
        assert len(tool.cli_command) > 0


def test_to_dict_is_serializable():
    for tool in REGISTRY:
        d = tool.to_dict()
        dumped = json.dumps(d)
        reloaded = json.loads(dumped)
        assert reloaded["name"] == tool.name


def test_to_json_global():
    result = to_json()
    parsed = json.loads(result)
    assert isinstance(parsed, list)
    assert len(parsed) == len(REGISTRY)


def test_get_tool_by_name():
    tool = get_tool("SeqNado")
    assert tool is not None
    assert tool.pypi_name == "seqnado"


def test_get_tool_by_pypi_name():
    tool = get_tool("plotnado")
    assert tool is not None
    assert tool.name == "PlotNado"


def test_get_tool_case_insensitive():
    assert get_tool("SEQNADO") == get_tool("seqnado")


def test_get_tool_missing_returns_none():
    assert get_tool("doesnotexist") is None


def test_python_compatible_known_cap():
    tool = get_tool("tabnado")
    assert tool is not None
    # python_compatible() returns a bool — just check it doesn't raise
    result = tool.python_compatible()
    assert isinstance(result, bool)


def test_bamnado_has_cli_command():
    bamnado = get_tool("bamnado")
    assert bamnado is not None
    assert bamnado.cli_command == "bamnado"


def test_regulonado_has_cli_command():
    regulonado = get_tool("regulonado")
    assert regulonado is not None
    assert regulonado.name == "ReguloNado"
    assert regulonado.cli_command == "regulonado"
    assert "bigWig" in regulonado.input_types


def test_library_only_tools_have_no_cli_command():
    mccnado = get_tool("mccnado")
    assert mccnado is not None and mccnado.cli_command is None


def test_input_output_types_are_lists():
    for tool in REGISTRY:
        assert isinstance(tool.input_types, list)
        assert isinstance(tool.output_types, list)
        assert len(tool.input_types) > 0
        assert len(tool.output_types) > 0
