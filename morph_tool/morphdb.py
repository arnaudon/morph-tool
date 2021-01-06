'''Provides two classes to work with neurondb files

The main use is to provide a neuronDB loader and a method to retrieve
information as a dataframe:

MorphologyDB('neurondb.xml').df
'''
import json
import logging
from pathlib import Path
from typing import Any, Callable, Dict, Iterable, List, Optional, Tuple

import pandas as pd
import xmltodict

L = logging.getLogger(__name__)

NEURONDB_XML = "neuronDB.xml"
NEURONDB_DAT = "neurondb.dat"
MTYPE_MSUBTYPE_SEPARATOR = ':'
EXTS = {".asc", ".ASC", ".h5", ".H5", ".swc", ".SWC"}


def find_morph(folder: Path, stem: str, ext: str) -> Optional[Path]:
    """Returns the path to a morphology in morphology_dir matching stem.

    If no morphology is found, returns None. Provide ext for speed up
    """
    if ext:
        return folder / (stem + ext)

    for ext in EXTS:
        path = folder / (stem + ext)
        if path.exists():
            return path
    return None


class MorphInfo:
    '''A class the contains information about a morphology.

    Its role is to abstract away the raw data.
    '''

    # Attributes found in the <repair></repair> block of the XML
    BOOLEAN_REPAIR_ATTRS = [
        'use_axon',
        'use_dendrites',
        'axon_repair',
        'dendrite_repair',
        'basal_dendrite_repair',
        'tuft_dendrite_repair',
        'oblique_dendrite_repair',
        'unravel',
        'use_for_stats',
    ]

    COLUMNS = (
        ["name", "mtype", "msubtype", "fullmtype", "layer", "label", "path"]
        + BOOLEAN_REPAIR_ATTRS
        + ["axon_inputs"]
    )

    def __init__(self, item: Dict[str, Any]):
        '''MorphInfo ctor.

        Args:
            item: A dictionary that represents the content of the XML file.
                  The only mandatory keys are [name, mtype, layer]
        '''
        self.name = item['name']
        self.mtype = item['mtype']
        self.layer = str(item['layer'])
        self.msubtype = item.get('msubtype') or ''
        self.fullmtype = MTYPE_MSUBTYPE_SEPARATOR.join(filter(None, [self.mtype, self.msubtype]))
        self.layer = str(item.get("layer", None))
        self.msubtype = item.get("msubtype") or ""

        def is_true(el: Optional[str]) -> bool:
            '''Parse a string representing a boolean repair flag and returns its boolean value

            Unless clearly stated as false, missing tags default to True

            - According to Eilif, an empty use_axon (corresponding to a null in the database)
              means that the axon is supposed to be used

            - dendrite_repair       defaults to True in BlueRepairSDK
            - basal_dendrite_repair defaults to True in BlueRepairSDK
            - unravel: well I guess we always want to do it
            '''
            return el in (None, '', 'true', 'True', )

        repair = item.get('repair', {})
        for attr in MorphInfo.BOOLEAN_REPAIR_ATTRS:
            setattr(self, attr, is_true(repair.get(attr)))

        # Case where there is a single <axoninput></axoninput> tag in the XML
        self.axon_inputs = repair.get('axon_sources', {}).get('axoninput', [])
        if not isinstance(self.axon_inputs, list):
            self.axon_inputs = [self.axon_inputs]

        # lineage information
        self.dendrite_donor = item.get('parent') or item.get('dendrite')
        self.axon_donor = item.get('parent') or item.get('axon')

    @property
    def data(self) -> Dict:
        '''Data that matter to generate the neurondb.xml'''
        return {
            'name': self.name,
            'mtype': self.mtype,
            'msubtype': self.msubtype,
            'layer': self.layer,
            'repair': {attr: getattr(self, attr)
                       for attr in MorphInfo.BOOLEAN_REPAIR_ATTRS + ['axon_inputs']}
        }

    @property
    def row(self) -> List:
        '''Flattened data structude ready to be used by a dataframe.'''
        return [getattr(self, attr) for attr in MorphInfo.COLUMNS]

    def __repr__(self):
        return f'MorphInfo(name={self.name}, mtype={self.mtype}, layer={self.layer})'


class MorphAPI(object):
class MorphologyDB(object):
    '''A MorphInfo container.

    It takes care of maintaining unicity of the MorphInfo element
    and methods to write neurondb to various format (xml, dat, csv)
    '''

    def __init__(self, root_path=None):
        """Builds a MorphologyDB from a morphology folder"""
        self.df = pd.DataFrame()
        if root_path:
            self.add_from_path(root_path)

    def _add_from_single_folder(self, path):
        """Load morphologies from a single folder with neurondb.xml."""

    def _add_from_folders(self, root_path):
        """Load several folders as a release with various neurondb.xml."""
        for folder in root_path.iterdir():
            if folder.is_dir():
                self.add_neurondb(folder / NEURONDB_XML)

    def add_from_path(self, root_path):
        """Load a morphology release."""
        root_path = Path(root_path).resolve()
        if (root_path / NEURONDB_XML).exists():
            self.add_neurondb(root_path / NEURONDB_XML)
        elif any((folder / NEURONDB_XML).exists() for folder in root_path.iterdir()):
            for folder in root_path.iterdir():
                if (folder / NEURONDB_XML).exists():
                    self.add_neurondb(folder / NEURONDB_XML)
        else:
            raise Exception(f"We cannot load morphologies from path {root_path}")

    def add_neurondb(
        self,
        neurondb: Path = None,
        label: str = None,
        ext=None,
        morphology_folder: Optional[Path] = None,
        morph_info_filter: Callable[[MorphInfo], bool] = None,
    ):
        """Builds a MorphologyDB from a path (several options available)

        Args:
            neurondb: path to a neurondb.xml/neurondb.dat or morphology folder
            label: a unique label to mark all morphologies coming from this neurondb
            morphology_folder: the location of the morphology files, if None it will default
                to the neurondb folder
            morph_info_filter: a filter function to be applied to each MorphInfo element
        """
        neurondb = Path(neurondb)
        loader_type = None
        if neurondb.is_dir():
            if (neurondb / NEURONDB_XML).exists():
                neurondb = neurondb / NEURONDB_XML
                loader_type = "xml"
            elif (neurondb / NEURONDB_DAT).exists():
                neurondb = neurondb / NEURONDB_DAT
                loader_type = "dat"
            else:
                # TODO: simple loading from a folder is some morphs exist
                raise Exception(f"We cannot load neurondb at {neurondb}")
        elif neurondb.suffix == ".xml":
            loader_type = "xml"
        elif neurondb.suffix == ".dat":
            loader_type = "dat"
        else:
            raise Exception(f"We cannot load neurondb at {neurondb}")

        if not morphology_folder:
            morphology_folder = neurondb.parent.resolve()

        if not label:
            label = morphology_folder.stem

        if not ext:
            # this is used to bypass the check of file existence later, which is slow
            last_split = "." + morphology_folder.stem.split("-")[-1]
            if last_split in EXTS:
                ext = last_split

        if loader_type == "xml":
            self._add_from_xml(
                neurondb, label, ext, morphology_folder, morph_info_filter
            )
        if loader_type == "dat":
            self._add_from_dat(
                neurondb, label, ext, morphology_folder, morph_info_filter
            )

    def _add_from_dat(
        self,
        neurondb: Path = None,
        label: str = None,
        ext=None,
        morphology_folder: Optional[Path] = None,
        morph_info_filter: Callable[[MorphInfo], bool] = None,
    ):
        """Builds a MorphologyDB from a neurondb.dat file

        TO IMPLEMENT
        """
        pass

    def _add_from_xml(
        self,
        neurondb: Path = None,
        label: str = None,
        ext=None,
        morphology_folder: Optional[Path] = None,
        morph_info_filter: Callable[[MorphInfo], bool] = None,
    ):
        """Builds a MorphologyDB from a neurondb.xml file

        Args:
            neurondb: path to a neurondb.xml file
            label: a unique label to mark all morphologies coming from this neurondb
            morphology_folder: the location of the morphology files, if None it will default
                to the neurondb folder
            morph_info_filter: a filter function to be applied to each MorphInfo element
        """
        with open(neurondb) as fd:
            content = fd.read()
            neurondb = xmltodict.parse(content)

        morphologies = neurondb["neurondb"]["listing"]["morphology"]
        for morph in morphologies:
            morph["label"] = label

        # Case where there is a single <morphology></morphology> tag in the XML
        if not isinstance(morphologies, list):
            morphologies = [morphologies]
        morphologies = filter(None, morphologies)
        morphologies = map(MorphInfo, morphologies)

        if morph_info_filter is not None:
            morphologies = filter(morph_info_filter, morphologies)

        morphologies = list(morphologies)
        for morph in morphologies:
            morph.path = find_morph(morphology_folder, morph.name, ext)

        dataframe = pd.DataFrame(
            [morph.row for morph in morphologies], columns=MorphInfo.COLUMNS
        )
        self.df = pd.concat([self.df, dataframe])

    @property
    def labels(self):
        """Return list of labels."""
        return self.df.label.unique()

    def by_label(self, label):
        """Return df with specific label."""
        return self.df.groupby("label").get_group(label)

    def add_morph(self, morph_info: MorphInfo) -> bool:
        """Add a morphology to the database."""
        self.df = self.df.append(morph_info.data).drop_duplicates.reset_index(drop=True)

    def add_morphs(self, morph_infos: list) -> bool:
        """Add a morphology to the database."""
        self.df = self.df.append(
            [morph_info.data for morph_info in morph_infos]
        ).drop_duplicates.reset_index(drop=True)

    def write(self,
              output_path: Path,
              filt: Dict[str, bool] = None):
        '''Write the neurondb file to XML, DAT or CSV format'''
        output_path = Path(output_path)
        ext = output_path.suffix.lower()[1:]
        if ext == "csv":
            self.df.to_csv(output_path, index=False)
        elif ext == "dat":
            self.df.to_csv(output_path, index=False, header=None, sep="\t")
        elif ext == "xml":
            with output_path.open("w") as fd:
                fd.write(xmltodict.unparse(self.data, pretty=True))
        else:
            raise ValueError(f'Unsupported neurondb extensions ({ext}).'
                             ' Should be one of: (xml,csv,dat)')
        return output_path

