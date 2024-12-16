from abc import ABC, abstractmethod


class FateWrapper(ABC):

    @abstractmethod
    def pipeline(self):
        pass

    def __contains__(self, item):
        "check if have attribute"
        return hasattr(self, item)

    def keys(self):
        """ return all attibute name, then the function dict() can be used"""
        return self.__dict__.keys()

    def __getitem__(self, key):
        "get attribute"
        return getattr(self, key)
