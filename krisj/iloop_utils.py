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
# Interface functions for the iLoop platform
# Kristian Jensen, October 2015
#

from __future__ import print_function

import requests
import json

try:
    from functools import lru_cache
except ImportError:
    def lru_cache(func):
        return func


class ILoopConnection(object):
    def __init__(self, root_url, auth_path, password, verify=True):
        self._root = root_url.rstrip("/")
        with open(auth_path) as infile:
            token = infile.read().strip()
        self._auth = (token, password)
        self._verify = bool(verify)
        if not self._verify:
            requests.packages.urllib3.disable_warnings()

    @staticmethod
    def _check_status(response):
        if response.status_code // 100 != 2:
            print(response.status_code)
            print(response.text)
        response.raise_for_status()

    def _make_get_request(self, path, data=None, params=None):
        r = requests.get(self._root+"/"+path.lstrip("/"), auth=self._auth, json=data, params=params, verify=self._verify)
        self._check_status(r)
        return r

    def _make_post_request(self, path, data=None):
        if data is None:
            data = {}
        r =  requests.post(self._root+"/"+path.lstrip("/"), auth=self._auth, json=data, verify=self._verify)
        self._check_status(r)
        return r

    def _make_delete_request(self, path, data=None, params=None):
        r = requests.delete(self._root+"/"+path.lstrip("/"), auth=self._auth, json=data, params=params, verify=self._verify)
        self._check_status(r)
        return r

    def _make_patch_request(self, path, data=None, params=None):
        if data is None:
            data = {}
        r = requests.patch(self._root+"/"+path.lstrip("/"), auth=self._auth, json=data, params=params, verify=self._verify)
        self._check_status(r)
        return r

    def _add_object(self, object_type, object_data):
        return self._make_post_request("/api/"+object_type, object_data)

    def _add_recipe_to_medium(self, medium_id, recipe):
        data = []
        for k, v in recipe.items():
            data.append({"compound": k, "concentration": v})
        self._make_post_request("/api/medium/"+str(medium_id)+"/contents", data)

    def add_medium(self, name, identifier, ph, description=None, composition=None):
        if composition is not None:  # Check that all compounds are in the DB
            for chebi_name in composition:
                result = self.list_compounds("chebi_name", chebi_name)
                if len(result) == 0:
                    raise ValueError(chebi_name + " could not be found in the database.")
        data = {"name": name, "identifier": identifier, "ph": ph, "description": description}
        r = self._add_object("medium", data)
        medium_id = self._get_id(r.json())
        if composition is not None:
            self._add_recipe_to_medium(medium_id, composition)

    def compound_to_chebi(self, compound_name):
        r = self._make_get_request("/api/chemical-entity/search", {"query": compound_name})
        data = r.json()
        if data[0]["chebi_name"] == compound_name:
            return data[0]["chebi_id"]
        else:
            raise ValueError(compound_name + " was not found.")

    def add_pool(self, identifier, project, alias, description=None, **kwargs):
        project = self._get_project_id(project)
        data = {
            "identifier": identifier,
            "project": project,
            "description": description,
            "alias": alias}
        data.update(kwargs)
        r = self._add_object("pool", data)
        return r.json()

    def add_organism(self, name, short_code, alternate_name=None):
        r = self._add_object(
            "organism", {"name": name, "short_code": short_code, "alternate_name": alternate_name})
        return r.json()

    @lru_cache(maxsize=None)
    def _get_organism_id(self, short_code):
        r = self._make_get_request("/api/organism", params=self._where("short_code", short_code))
        try:
            organism = r.json()[0]
        except IndexError:
            raise ValueError("Organism not found")
        return self._get_id(organism)

    @lru_cache(maxsize=None)
    def _get_project_id(self, project_code=None, project_name=None):
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
        return self._get_id(project)

    def add_strain(self, identifier, alias, organism, project, description=None, **kwargs):
        organism = self._get_organism_id(organism)
        project = self._get_project_id(project)
        data = {"identifier": identifier,
                "alias": alias,
                "organism": organism,
                "description": description,
                "project": project}
        data.update(kwargs)
        r = self._add_object("strain", data)
        return r.json()

    def add_plate(self, barcode, project, type, contents=None):
        project = self._get_project_id(project)
        data = {"barcode": barcode,
                "project": project,
                "type": type}
        r = self._add_object("plate", data)
        plate = r.json()
        if contents is not None:
            self._make_post_request(plate["$uri"]+"/contents", contents)
        return plate

    def add_experiment(self, identifier, project, type, title=None, device=None, **kwargs):
        project = self._get_project_id(project)
        data = {"identifier": identifier,
                "project": project,
                "type": type,
                "device": device,
                "title": title}
        data.update(kwargs)
        r = self._add_object("experiment", data)
        return r.json()

    def add_sample(self, experiment, strain=None, medium=None,
                   plate=None, position=None, scalars=None, series=None):
        data = [{"strain": strain,
                "medium": medium,
                "plate": plate,
                "position": position,
                "scalars": scalars,
                "series": series
        }]
        return self.add_samples(experiment, data)

    def add_samples(self, experiment, data):
        """Data is a list of dicts. Possible data fields are 'strain', 'medium', 'plate', 'position',
        'scalars', 'series'.
        """
        uri = experiment["$uri"]+"/samples"
        r = self._make_post_request(uri, data)
        return r.json()

    @staticmethod
    def _where(field, query, type="is"):
        try:
            query = {"$ref": query["$uri"]}
        except (TypeError, KeyError):
            pass
        query_types = ["startswith", "endswith", "ne"]
        if type == "is":
            return {"where": json.dumps({field: query})}
        elif type in query_types:
            return {"where": json.dumps({field: {"$"+type: query}})}
        else:
            raise ValueError("Invalid type. Choose from "+str(query_types))

    @staticmethod
    def _get_id(obj):
        return int(obj["$uri"].split("/")[-1])

    def _list_object(self, object_type, *args):
        if args:
            params = self._where(*args)
        else:
            params = None
        r = self._make_get_request("/api/"+object_type, params=params)
        return r.json()

    def list_media(self, *args):
        return self._list_object("medium", *args)

    def list_pools(self, *args):
        return self._list_object("pool", *args)

    def list_strains(self, *args):
        return self._list_object("strain", *args)

    def list_compounds(self, *args):
        return self._list_object("chemical-entity", *args)

    def list_plates(self, *args):
        return self._list_object("plate", *args)

    def _get_record_by(self, table, field, query, all=False):
        records = self._list_object(table, field, query)
        if all:
            return records
        else:
            if len(records) > 0:
                return records[0]
            else:
                raise ValueError(table.capitalize()+" '"+query+"' could not be found.")

    def get_organism(self, query, field="short_code", **kwargs):
        return self._get_record_by("organism", field, query, **kwargs)

    def get_project(self, query, field="code", **kwargs):
        return self._get_record_by("project", field, query, **kwargs)

    @lru_cache(maxsize=None)
    def get_compound(self, query, field="chebi_name", **kwargs):
        return self._get_record_by("chemical-entity", field, query, **kwargs)

    def get_strain(self, query, field="alias", **kwargs):
        return self._get_record_by("strain", field, query, **kwargs)

    @property
    def plate_types(self):
        r = self._make_get_request("/api/plate/types")
        return [p["type"] for p in r.json()]
