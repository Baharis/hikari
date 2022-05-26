"""This file contains tools for making property maps visualised on sphere"""


class AngularPropertyExplorerFactory:
    """Factory method for angular prop. explorers based on realpython thread."""
    def __init__(self):
        self._explorers = {}

    def register_explorer(self, prop, explorer):
        self._explorers[prop] = explorer

    def create(self, prop, **kwargs):
        explorer = self._explorers.get(prop)
        if not explorer:
            raise ValueError(f'Explorer for {prop} has not been registered!')
        return explorer(**kwargs)
