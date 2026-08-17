"""配布データから必要なぶんだけ取り出す."""
import os

import requests

BASE = "."
DATA = "data"
TRACKS = "tracks"


def setup(base, data="data", tracks="tracks"):
    """取得先と作業ディレクトリを決める. base は URL でも手元のパスでもよい."""
    global BASE, DATA, TRACKS
    BASE, DATA, TRACKS = base, data, tracks
    for directory in (DATA, TRACKS):
        os.makedirs(directory, exist_ok=True)


def fetch(name, first=None, last=None):
    """name のバイト列を返す. first と last を渡すとその範囲だけ.

    範囲取得の結果は DATA に控えるので, 同じ範囲を二度落とさない.
    """
    cached = None if first is None else f"{DATA}/.range_{name}_{first}_{last}"
    if cached and os.path.exists(cached):
        return open(cached, "rb").read()

    if BASE.startswith("http"):
        headers = {"Range": f"bytes={first}-{last - 1}"} if first is not None else {}
        response = requests.get(f"{BASE}/{name}", headers=headers, timeout=600)
        response.raise_for_status()
        data = response.content
    else:
        with open(f"{BASE}/{name}", "rb") as fh:
            if first is None:
                data = fh.read()
            else:
                fh.seek(first)
                data = fh.read(last - first)

    if cached:
        open(cached, "wb").write(data)
    return data


def download(name):
    """ファイル全体を DATA に置き, そのパスを返す. すでにあれば何もしない."""
    path = f"{DATA}/{name}"
    if not os.path.exists(path):
        print(f"downloading {name}")
        with open(path, "wb") as fh:
            fh.write(fetch(name))
    return path
