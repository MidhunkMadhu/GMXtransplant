"""Field definitions derived from the pipeline schema, without GUI dependencies."""
from dataclasses import asdict, is_dataclass
from functools import lru_cache
from typing import Any, Dict, Optional, Union, get_args, get_origin, get_type_hints
import config
import types
from charmprot import CharmProtSpec
from addbinder import AddBinderSpec


@lru_cache(None)
def hints(cls):
    return get_type_hints(cls)


def alternatives(annotation):
    return get_args(annotation) if get_origin(annotation) in (Union, types.UnionType) else (annotation,)


def concrete(annotation):
    return next((t for t in alternatives(annotation) if t is not type(None)), str)


def annotation_for(path):
    if path[:2] == ('minimization', 'resources') and len(path) == 3:
        return int if path[-1] == 'cpus_per_task' else str
    annotation = config.Config
    for index, key in enumerate(path):
        annotation = concrete(annotation)
        if index == 0 and key == 'charmprot':
            annotation = CharmProtSpec
        elif index == 0 and key == 'addbinder':
            annotation = AddBinderSpec
        elif index == 0 and key == 'paths':
            annotation = Dict[str, Optional[str]]
        elif is_dataclass(annotation):
            annotation = hints(annotation).get(key, Any)
        elif get_origin(annotation) is dict:
            args = get_args(annotation)
            annotation = args[1] if args else Any
        elif get_origin(annotation) is list:
            args = get_args(annotation)
            annotation = args[0] if args else Any
        else:
            return Any
    return annotation


def initial_value(annotation, path=()):
    annotation = concrete(annotation)
    if is_dataclass(annotation):
        return asdict(annotation())
    origin = get_origin(annotation)
    if origin is list:
        if path[-1:] == ('box_dimensions',):
            return [0.0, 0.0, 0.0, 90.0, 90.0, 90.0]
        return []
    if origin is dict or annotation is dict:
        if path[:2] == ('cholesterol', 'composition') and len(path) == 4 and path[2] == 'lipid_targets':
            return {'upper': 0, 'lower': 0}
        return {}
    return {bool: False, int: 0, float: 0.0}.get(annotation, '')


def residue_mapping(path):
    return path in {('name_restoration', 'pdb_to_full_resname'),
                    ('name_restoration', 'topology_to_pdb_resname'),
                    ('cholesterol', 'composition', 'lipid_targets'),
                    ('topology', 'moleculetype_overrides'),
                    ('charge', 'charge_table_overrides'),
                    ('minimization', 'restraint_residue_classes')}
