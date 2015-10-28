# Copyright 2015 Novo Nordisk Foundation Center for Biosustainability, DTU.

# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at

# http://www.apache.org/licenses/LICENSE-2.0

# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

#
# Various utility functions
# Kristian Jensen, March 2015
#

from __future__ import print_function

import requests
import json


class ILoopConnection(object):
    def __init__(self, root_url, auth_path):
        self.root = root_url.rstrip("/")
        with open(auth_path) as infile:
            token = infile.read().strip()
        self.auth = (token, "x-oauth-basic")

    def _make_get_request(self, path, data=None, params=None):
        if data is None:
            data = {}
        return requests.get(self.root+"/"+path.lstrip("/"), auth=self.auth, json=data, params=params)

    def _make_post_request(self, path, data):
        return requests.post(self.root+"/"+path.lstrip("/"), auth=self.auth, json=data)

    def _add_object(self, object_type, object_data):
        return self._make_post_request("/api/"+object_type, object_data)

    def _add_recipe_to_medium(self, medium_id, recipe):
        data = []
        for k, v in recipe.items():
            data.append({"compound": k, "concentration": v})
        self._make_post_request("/medium/"+str(medium_id)+"/contents", data)

    def add_medium(self, name, identifier, ph, description=None, composition=None):
        data = {"name": name, "identifier": identifier, "ph": ph, "description": description}
        r = self._add_object("medium", data)
        medium_id = r.json()["$uri"].split("/")[-1]
        if composition is not None:
            self._add_recipe_to_medium(medium_id, composition)

    def compound_to_id(self, compound_name):
        r = self._make_get_request("/api/chemical-entity/search", {"query": compound_name})
        data = r.json()
        if data[0]["chebi_name"] == compound_name:
            return data[0]["chebi_id"]
        else:
            raise ValueError(compound_name + " was not found.")

    def add_pool(self, identifier, project, description=None):
        r = self._add_object("pool", {
            "identifier": identifier, "project": project, "description": description})
        return r.json()

    def add_organism(self, name, short_code, alternate_name=None):
        r = self._add_object(
            "organism", {"name": name, "short_code": short_code, "alternate_name": alternate_name})
        return r.json()

    def _get_organism_uri(self, short_code):
        r = self._make_get_request("/api/organism", params=self._where("short_code", short_code))
        try:
            organism = r.json()[0]
        except IndexError:
            raise ValueError("Organism not found")
        return organism["$uri"]

    def _get_project_uri(self, project_code=None, project_name=None):
        if project_code is not None:
            q = project_code
            field = "code"
        elif project_name is not None:
            q = project_name
            field = "name"
        else:
            raise ValueError("You must supply a project code or a project name")
        r = self._make_get_request("/api/project", params=self._where(field, q))
        try:
            project = r.json()[0]
        except IndexError:
            raise ValueError("Project not found")
        return project["$uri"]

    def add_strain(self, identifier, alias, organism, description=None):
        organism = self._get_organism_uri(organism)
        r = self._add_object("strain", {
                             "identifier": identifier,
                             "alias": alias,
                             "organism": organism,
                             "description": description})
        return r.json()

    @staticmethod
    def _where(field, query):
        return {"where": json.dumps({field: query})}
