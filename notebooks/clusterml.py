#!/usr/bin/env python

# %%
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import gzip
import json
import qml


class MockXYZ(object):
    def __init__(self, lines):
        self._lines = lines

    def readlines(self):
        return self._lines


# %%
# read source data, remove diradicals
df = pd.read_csv(
    "../data/qm9-shape-properties",
    sep=" ",
    names="dbid dipole polarizability homo lumo gap extent zpve internalzero internalfinite enthalpy freeenergy heatcap volume area".split(),
)
with open("../data/uncharacterized.txt") as fh:
    lines = fh.readlines()[9:-1]
    blackout = [int(_.strip().split()[0]) for _ in lines]
df = df[~df.dbid.isin(blackout)]

# %%
# annotate with shapes
def classify(volume, area):
    def spheroid(xs, ratio):
        e = 1 - 1 / ratio ** 2
        return (
            2
            * np.pi
            * (3 * xs / (4 * np.pi * ratio)) ** (2 / 3)
            * (1 + ratio / np.sqrt(e) * np.arcsin(np.sqrt(e)))
        )

    def cuboid(xs, px, py):
        return 2 * (px + py + px * py) * (xs / (px * py)) ** (2 / 3)

    def disk(xs, ratio):
        return 2 * np.pi * (xs * ratio / np.pi) ** (2 / 3) * (1 + 1 / ratio)

    def vertical(xs, limit):
        return [0] + [limit] * (len(xs) - 1)

    def eldisk(xs, px, py):
        h = (1 - px) ** 2 / (1 + px) ** 2
        q = 2 * np.pi * px + py * np.pi * (1 + px) * (
            1 + 3 * h / (10 + np.sqrt(4 - 3 * h))
        )
        return (xs / (np.pi * px * py)) ** (2 / 3) * q

    def sphere(xs):
        return (6 * xs) ** (2 / 3) * (np.pi) ** (1 / 3)

    def cube(xs):
        return 6 * (xs) ** (2 / 3)

    if volume < 10 and area < 10:
        return "linear"

    borders = np.array(
        [
            eldisk(volume, 25, 0.2),
            eldisk(volume, 8, 0.2),
            eldisk(volume, 1.5, 0.2),
            disk(volume, 5),
            disk(volume, 3),
            cuboid(volume, 1, 5),
            cuboid(volume, 1, 2),
            cube(volume),
            spheroid(volume, 3),
            spheroid(volume, 1.5),
            sphere(volume),
        ]
    )
    labels = np.array(
        [
            "planar",
            "eldisk",
            "eldisk",
            "disk",
            "disk",
            "cuboid",
            "cuboid",
            "cube",
            "spheroid",
            "spheroid",
            "sphere",
        ]
    )

    return labels[np.argmin(np.abs(borders - area))]


df["shape"] = df.apply(lambda row: classify(row.volume, row.area), axis=1)

# %%
# read geometries
with gzip.open("../../rank-only/src/.QM9.cache.gz", "r") as fh:
    data = fh.read()
    cache = json.loads(data.decode("ascii"))
FILENAMES = np.array(cache["filenames"])
GEOMETRIES = np.array(cache["geometries"])
DBIDS = [int(_.split("_")[1].split(".")[0].lstrip("0")) for _ in FILENAMES]

# %%
def get_representation(repname, dbid):
    xyz = MockXYZ(GEOMETRIES[DBIDS.index(dbid)].split("\n"))
    if repname == "CM":
        mol = qml.Compound(xyz=xyz)
        mol.generate_coulomb_matrix(size=34, sorting="row-norm")
        return mol.representation
    if repname == "BoB":
        mol = qml.Compound(xyz=xyz)
        mol.generate_bob(size=23, asize={"C": 6, "H": 20, "N": 2, "O": 2, "F": 2})
        return mol.representation


# %%
def get_dataset(database, ids, propertyname, repname):
    s = pd.merge(pd.DataFrame({"dbid": ids}), database, how="left").sort_values(
        by="dbid"
    )
    X = np.array([get_representation(repname, _) for _ in s.dbid.values])
    Y = s[propertyname].values
    return X, Y


def do_split(database, N, propertyname, repname):
    allids = database.dbid.values
    np.random.shuffle(allids)
    trainids, testids = allids[:N], allids[N:]
    trainX, trainY = get_dataset(database, trainids, propertyname, repname)
    testX, testY = get_dataset(database, testids, propertyname, repname)
    return trainX, trainY, testX, testY


def learn_and_predict(database, N, repname, sigma, lval, propname):
    trainX, trainY, testX, testY = do_split(database, N, propname, repname)

    K = qml.kernels.gaussian_kernel(trainX, trainX, sigma)
    K[np.diag_indices_from(K)] += lval
    alphas = qml.math.cho_solve(K, trainY)

    K = qml.kernels.gaussian_kernel(trainX, testX, sigma)
    Ys = np.dot(K.transpose(), alphas)

    return np.abs(Ys - testY).mean()


# %%
learn_and_predict(df.query("shape=='planar'"), 100, "CM", 2, 1e-7, "lumo")

# %%
