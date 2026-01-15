from typing import Callable, Any
from os.path import exists
from os import makedirs
import json


class JsonManager:
    def __init__(
        self,
        savepath: str,
        complete_override: bool = False,
        override_files: list[str] | None = None,
    ) -> None:
        if complete_override:
            print("Complete override set to True. Proceed?")
            if input() != "y":
                raise Exception("Program terminated.")
        self.savepath = savepath
        self.override_files = override_files
        self.complete_override = complete_override
        self.file_to_computation: dict[str, Callable[[], dict[str, Any]]] = {}
        self.__cached_data: dict[str, dict[str, Any]] = {}

    def check_folder(self) -> None:
        if not exists(self.savepath):
            makedirs(self.savepath)

    def _check_json(self, json_file: str | None = None) -> bool:
        self.check_folder()
        return exists(self.savepath + f"/{json_file}.json")

    def get_data_from_json(self, filename: str, key: str) -> Any:
        """
        Retrieves from the json file 'filename.json' the
        value associated to the input key.
        """
        if filename not in self.__cached_data:
            self.__cached_data[filename] = {}
        if key in self.__cached_data[filename]:
            return self.__cached_data[filename][key]
        if (
            not self._check_json(filename)
            or (self.override_files and filename in self.override_files)
            or self.complete_override
        ):
            data = self.file_to_computation[filename]()
            with open(f"{self.savepath}/{filename}.json", "w") as f:
                json.dump(data, f)
            return data[key]
        with open(f"{self.savepath}/{filename}.json", "r") as f:
            data = json.load(f)
            self.__cached_data[filename].update(data)
            return data[key]
