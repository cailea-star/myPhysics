import json
import sys
from pathlib import Path

from sklearn.feature_extraction.text import TfidfVectorizer
from sklearn.metrics.pairwise import cosine_similarity


ROOT_PATH = Path(__file__).resolve().parent.parent
TAGS_PATH = ROOT_PATH / "vocab" / "tags.json"


def load_tags(path: Path = TAGS_PATH) -> list[dict]:
    return json.loads(path.read_text(encoding="utf-8"))


def tag_normalize(tag: dict) -> str:
    parts = [tag["tag"], *tag.get("aliases", []), tag.get("definition", "")]
    return " ".join(parts).replace("_", " ").replace("-", " ").casefold()


def tag_score(query_vector, tag_vector) -> float:
    return float(cosine_similarity(query_vector, tag_vector)[0, 0])


def search_similar_tags(query: str, tags: list[dict]) -> list[tuple[float, dict]]:
    documents = [query.casefold(), *(tag_normalize(tag) for tag in tags)]
    vectors = TfidfVectorizer(analyzer="char_wb", ngram_range=(2, 4)).fit_transform(documents)
    matches = [(tag_score(vectors[0], vectors[i]), tag) for i, tag in enumerate(tags, 1)]
    return sorted(matches, key=lambda match: (-match[0], match[1]["tag"]))


def tag_print(score: float, tag: dict) -> None:
    print(
        f"{score:.3f}\t{tag['tag']}\t{', '.join(tag.get('types', []))}\t"
        f"{', '.join(tag.get('aliases', []))}\t{tag.get('definition', '')}"
    )


def main(query: str, number: int) -> int:
    if number < 1:
        raise SystemExit("number must be a positive integer")
    for score, tag in search_similar_tags(query, load_tags())[:number]:
        tag_print(score, tag)
    return 0


if __name__ == "__main__":
    if len(sys.argv) != 3:
        raise SystemExit("usage: python scripts/search_similar_tags.py QUERY NUMBER")
    raise SystemExit(main(sys.argv[1], int(sys.argv[2])))
