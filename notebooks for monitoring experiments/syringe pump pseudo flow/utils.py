from dataclasses import dataclass
import os
from networkx import NetworkXNoPath, shortest_path_length

from AnalyticalLabware import SpinsolveNMR, SpinsolveNMRSpectrum
from chempiler import Chempiler
from chempiler.chempiler import ChempilerGraph


class DiscoveryPlatform:
    def __init__(self, nmr: SpinsolveNMR, c: Chempiler):
        self.nmr = nmr
        self.c = c

    # NMR related

    def test_nmr(self):
        self.nmr.sample = "test"
        self.nmr.get_spectrum(("1D PROTON", {"Scan": "QuickScan"}))

    def shim_on_sample(self, ref: float):
        self.nmr.sample = "shim"
        self.nmr.get_spectrum(
            (
                "SHIM 1H SAMPLE",
                {
                    "Mode": "Manual",
                    "manualStart": -25,
                    "manualEnd": 35,
                    "SampleReference": ref,
                    "Shim": "QuickShim1",
                },
            )
        )

    def get_wet_sup(self):
        self.nmr.get_spectrum(
            (
                "1D WET SUP",
                {
                    "Mode": "Manual",
                    "satFrequency1": 1.6,
                    "CorrectionFactor": 1,
                    "Dummy": 0,
                    "Number": 2,
                    "AcquisitionTime": 3.2,
                    "RepetitionTime": 4,
                },
            )
        )

    def get_fluorine(self):
        self.nmr.get_spectrum(("1D FLUORINE+", {"Number": "128"}))

    # Chempiler related

    def take_sample(self, pump: str):
        self.c.move("sample", pump, 1)

    def return_sample(self, pump: str):
        self.c.move(pump, "sample", 1)

    def flowcell(self, pump1: str, pump2: str, volume: float = 1):
        self.c.move(
            pump1,
            pump2,
            volume,
            speed=20,
            through_nodes=["in_splitter", "nmr", "out_splitter"],
        )

    def purge_flowcell(self, vessel: str | None = None, repeats: int = 1):
        vessel = vessel or "w23"
        for _ in range(repeats):
            self.c.move(
                "a22",
                vessel,
                5,
                initial_pump_speed=150,
                mid_pump_speed=20,
                end_pump_speed=150,
                through_nodes=["out_splitter", "nmr", "in_splitter"],
            )

    def clean_flowcell(self, vessel: str, repeats: int = 1, shim: bool = False):
        for _ in range(repeats):
            self.c.move(vessel, "p33", 10 if shim else 5, speed=100)
            self.flowcell("p33", "p32", 5)
            if shim:
                self.shim_on_sample(1.8)
                self.c.move("p33", "w23", 5, speed=150)
            self.c.move("p32", "w22", 5, speed=150)
            self.purge_flowcell()

    def sample_flask(self, flask: str, exp_code: str | None = None):
        self.nmr.sample = flask if exp_code is None else exp_code + "_" + flask
        self.c.move(flask, "p33", 10)
        self.flowcell("p33", "p32")
        # self.get_wet_sup()
        self.nmr.get_spectrum()
        self.flowcell("p32", "p33")
        self.c.move("p33", flask, 10)
        self.purge_flowcell(flask)
        self.clean_flowcell("MeCN", 2)

    def get_latest_spectrum(self) -> SpinsolveNMRSpectrum:
        spec = self.nmr.spectrum
        spec.default_processing()
        spec.autophase()
        spec.correct_baseline()
        spec.reference_spectrum(new_position=1.8, reference="highest")
        spec.normalise()

        return spec


def get_nearest_node(graph: ChempilerGraph, src: str, target_class: str):
    short_p_len = 100
    nearest_node = None
    target_nodes = [
        node for node in graph.nodes if graph.nodes[node]["class"] == target_class
    ]
    for node in target_nodes:
        try:
            p_len = shortest_path_length(graph, src, node)
        except NetworkXNoPath:
            try:
                p_len = shortest_path_length(graph, node, src)
            except NetworkXNoPath:
                continue
        if p_len < short_p_len:
            short_p_len = p_len
            nearest_node = node
    return nearest_node


def add_spectra(
    sm1: SpinsolveNMRSpectrum, sm2: SpinsolveNMRSpectrum
) -> SpinsolveNMRSpectrum:
    combined_spectrum = SpinsolveNMRSpectrum()
    combined_spectrum.x = sm1.x
    combined_spectrum.y = sm1.y + sm2.y
    return combined_spectrum


@dataclass
class Reaction:
    reagents: list[str]
    reactivity: float


def get_folder_paths(parent_folder):
    folder_paths = []
    for root, dirs, _ in os.walk(parent_folder):
        for directory in dirs:
            if "processed" not in directory and "Enhanced" not in directory:
                folder_paths.append(os.path.join(root, directory))
    return folder_paths
