import pytest
import cfe

import os
import yaml

import pandas as pd


class TestBackend:
    """abstract, cannot be instantiated"""
    pass


class TestDefinition():

    def setup_method(self):
        # definition in python function backend
        definition_yaml_filename = f"{os.path.dirname(__file__)}/../../cfe/method/definition/cf_paga.yml"
        with open(definition_yaml_filename, 'r') as file:
            definition_raw = yaml.safe_load(file)
        self.definition = cfe.method.Definition(definition_raw)

    def test_magic_methods(self):
        definition = self.definition

        definition["run"]

        # test __contains__
        assert "method" in definition, f"method should in definition"

        # test __getitem__
        assert definition["method"] == definition.method, f"definition['method'] should be the same as definition.method"

        # test keys
        definition_dict = dict(definition)
        attribute_name_list = ["method", "wrapper", "container", "package", "manuscript", "parameters"]
        assert set(attribute_name_list).issubset(set(definition_dict.keys())), f"{attribute_name_list} should be the keys of the dict: {definition_dict}"

    def test_get_inputs_df(self):
        inputs_df = self.definition.get_inputs_df()
        assert isinstance(inputs_df, pd.DataFrame)
        assert inputs_df.columns.tolist() == ["input_id", "required", "type"]
        assert inputs_df["type"].isin(["expression", "parameter", "prior_information"]).all()

    def test_get_parameters(self):
        parameters = self.definition.get_parameters()
        assert isinstance(parameters, dict)


if __name__ == "__main__":
    pytest.main(["-v", __file__])
