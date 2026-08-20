## Evo 2 Sparse Autoencoder Handson

ゲノム基盤モデル Evo 2 の内部表現を, スパースオートエンコーダ (SAE) で解釈可能な特徴に分解して読む
ハンズオンです. 大腸菌ゲノムを題材に, rRNA・tRNA・遺伝子間・ORF に対応する特徴を 32768 個の中から
探し, ゲノムブラウザとタンパク質の立体構造の上に読み出します. Evo 2 論文 (Brixi et al. *Nature* 2026)
Fig. 4 の大腸菌の解析をたどる構成です.

ノートブックは部ごとに分かれていて, 今後も部を追加していきます. 実行結果は
[GitHub Pages](https://genaibio.github.io/evo2-sae-handson/) で読めます.

| | | |
|---|---|---|
| [第 1 部 ゲノムトラックと立体構造への読み出し](https://genaibio.github.io/evo2-sae-handson/part1.html) | [evo2_sae_part1.ipynb](notebooks/evo2_sae_part1.ipynb) 7 節 | [![Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/GenAIBio/evo2-sae-handson/blob/main/notebooks/evo2_sae_part1.ipynb) |
| [第 2 部 アノテーションによる特徴の選択](https://genaibio.github.io/evo2-sae-handson/part2.html) | [evo2_sae_part2.ipynb](notebooks/evo2_sae_part2.ipynb) 3 節 | [![Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/GenAIBio/evo2-sae-handson/blob/main/notebooks/evo2_sae_part2.ipynb) |
| [第 3 部 DSSP・配列書き換え・プロファージ間の比較](https://genaibio.github.io/evo2-sae-handson/part3.html) | [evo2_sae_part3.ipynb](notebooks/evo2_sae_part3.ipynb) 4 節 | [![Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/GenAIBio/evo2-sae-handson/blob/main/notebooks/evo2_sae_part3.ipynb) |
| [付録 左文脈・エンコード・活性行列](https://genaibio.github.io/evo2-sae-handson/appendix.html) | [evo2_sae_appendix.ipynb](notebooks/evo2_sae_appendix.ipynb) 4 節 | [![Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/GenAIBio/evo2-sae-handson/blob/main/notebooks/evo2_sae_appendix.ipynb) |

手元で動かす場合は `jupyter lab notebooks/evo2_sae_part1.ipynb` です. 依存パッケージは最初の
セルが入れます. データは Hugging Face から取得し, `data/` に残ります.

| 場所 | |
|---|---|
| `notebooks/` | ノートブック |
| `evo2_sae/` | 取得・ゲノムブラウザ・立体構造の描画. ノートブックが import する |
| `dataset/` | 配布データの説明 (`README.md`) とチェックサム |
| `docs/` | GitHub Pages 用に書き出した HTML |

### データ

対象は大腸菌 *Escherichia coli* K-12 MG1655 (NC_000913.3) です. Evo 2 (7B) の `blocks.26` から
取り出した埋め込みと, Goodfire の SAE (`Goodfire/Evo-2-Layer-26-Mixed`) を通した全ゲノムの活性
行列を [Hugging Face](https://huggingface.co/datasets/suzuki-2001/evo2-sae-handson) に置いて
あります. 行列には 1 塩基あたり上位 128 個を保存してあり, HTTP Range で必要なバイト範囲だけを
取り出せます. 同じものを [Zenodo](https://doi.org/10.5281/zenodo.21856795) にも
アーカイブしてあります.

各ファイルの内容と, 配列・構造の取得元とライセンスは
[dataset/README.md](dataset/README.md) にあります.
