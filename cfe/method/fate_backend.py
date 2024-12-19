from abc import ABC, abstractmethod
import pandas as pd

from .._logging import logger


# Backend: abstract class, used for subsequent specific implementation such as "DockerBackend" class and "FunctionBackend" class
class Backend(ABC):
    def __init__(self, method_name):
        pass

    @abstractmethod
    def load_backend(self):
        pass

    @abstractmethod
    def run(self):
        pass

    @abstractmethod
    def _load_definition(self):
        """
        Load definition yml file during load_backend
        """
        pass

    def _extract_inputs(self, fdata, inputs_df):
        """
        ref: PyDynverse/pydynverse/wrap/method_extract_args.py _method_extract_inputs
        """
        # logger.debug("FateMethod _extract_inputs")

        # extract model input expression matrix
        input_ids = inputs_df["input_id"][inputs_df["type"] == "expression"].tolist()
        inputs = {}
        for expression_id in input_ids:
            # inputs[expression_id] = get_expression(dataset, expression_id)
            inputs[expression_id] = fdata.layers[expression_id]
        # main expression matrix, for example, Component1 and Slingshot need "expression", while monocle_ddrtree need "counts"
        inputs["expression_id"] = input_ids[0]
        # add cell and gene ids
        inputs["cell_ids"] = fdata.obs.index.tolist()
        inputs["feature_ids"] = fdata.var.index.tolist()
        return inputs

    def _extract_priors(self, fdata, inputs_df):
        """
        ref: PyDynverse/pydynverse/wrap/method_extract_args.py _method_extract_priors
        """
        # logger.debug("FateMethod _extract_priors")

        # extract prior information from
        priors = fdata.prior_information
        priors_key_list = priors.keys()
        priors_key_set = set(priors_key_list)

        # check required priors
        required_prior_ids = inputs_df["input_id"][inputs_df["required"] & (inputs_df["type"] == "prior_information")].tolist()
        required_prior_ids_set = set(required_prior_ids)
        if not (required_prior_ids_set <= priors_key_set):
            # all required priors are needed, if not , raise error
            missing_priors = required_prior_ids_set - priors_key_set
            msg = f"""
                ! Prior information {','.join(missing_priors)} is missing from dataset {fdata.id} but is required by the method. \n
                -> If known, you can add this prior information using fadata.add_prior_information({' ,'.join([str(i)+' = <prior>' for i in missing_priors])}). \n
                -> Otherwise, this method cannot be used.
            """
            raise Exception(msg)
        required_prior = {k: priors[k] for k in required_prior_ids}

        # check optional priors
        optional_prior_ids = inputs_df["input_id"][(~inputs_df["required"]) & (inputs_df["type"] == "prior_information")].tolist()
        optional_prior_ids_set = set(optional_prior_ids)
        if not (optional_prior_ids_set <= priors_key_set):
            # all required priors are not needed, enven if not all are provided, only warning
            missing_priors = list(optional_prior_ids_set - priors_key_set)
            msg = f"""
                Prior information {','.join(missing_priors)} is optional, but missing from dataset {fdata.id}. \n
                Will not give this prior to method.
            """
            logger.warning(msg)
        optional_prior = {k: priors[k] for k in list(optional_prior_ids_set & priors_key_set)}  # remove irrelevant keys

        priors = required_prior | optional_prior

        return priors


class Definition():

    def __init__(self, definition_raw: dict):
        self.method = definition_raw["method"]
        self.wrapper = definition_raw["wrapper"]
        self.container = definition_raw["container"]
        self.package = definition_raw["package"] if "package" in definition_raw else None
        self.manuscript = definition_raw["manuscript"] if "manuscript" in definition_raw else None
        self.parameters = pd.DataFrame(definition_raw["parameters"]).set_index("id")
        self.run = {}

        # inputs
        inputs = self.wrapper["input_required"]
        inputs = inputs.copy() if type(inputs) == list else [inputs]
        if "input_optional" in self.wrapper:
            input_optional = self.wrapper["input_optional"]
            inputs += input_optional if type(input_optional) == list else [input_optional]

        # extra input, including data and parameters
        params = self.parameters.index.tolist()
        input_id_list = inputs + params
        required_list = [i in self.wrapper["input_required"] for i in input_id_list]
        type_list = []
        for input_id in input_id_list:
            # type column
            if input_id in ["counts", "expression", "expression_future"]:
                type_list.append("expression")
            elif input_id in params:
                type_list.append("parameter")
            else:
                type_list.append("prior_information")
        inputs_df = pd.DataFrame({"input_id": input_id_list, "required": required_list, "type": type_list})
        self.wrapper["inputs"] = inputs_df

    def get_inputs_df(self):
        return self.wrapper["inputs"]

    def get_parameters(self):
        return dict(self.parameters["default"])

    def add_function_wrapper(self, return_function):
        if not return_function:
            # 直接返回字典格式
            return self
        # else:
        #     # 返回函数格式，等待默认参数进一步设置
        #     defaults = get_default_parameters(definition)  # 获取代码函数中的参数默认值

        #     def param_overrider_fun(**kwargs):
        #         # 参数覆盖, 接受代码函数中的参数默认值传入参数并覆盖definition.yml
        #         new_defaults = kwargs
        #         param_names = list(definition["parameters"].index)
        #         for param_name, v in new_defaults.items():
        #             if param_name in param_names:
        #                 definition["parameters"].loc[param_name, "default"] = v
        #             else:
        #                 # 该参数不definition.yml文件里定义
        #                 logger.error(f"Unknown parameter: {param_name}")
        #         return definition

        #     param_overrider_fun.__kwdefaults__ = defaults

        #     return param_overrider_fun
    def __contains__(self, item):
        "check if have attribute"
        return hasattr(self, item)

    def keys(self):
        """ return all attibute name, then the function dict() can be used"""
        return self.__dict__.keys()

    def __getitem__(self, key):
        "get attribute"
        return getattr(self, key)

    def __setitem__(self, key, value):
        "set attribute"
        setattr(self, key, value)
