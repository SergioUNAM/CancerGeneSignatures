from __future__ import annotations

import os
import time
from typing import Dict, List, Any, Optional, Tuple

import requests


class GoogleNLPClient:
    """Minimal client for Google Cloud Natural Language API (v1) using API key.

    Notes
    - Use API key via constructor, env GOOGLE_NLP_API_KEY or streamlit secrets.
    - Endpoints used: analyzeEntitySentiment, analyzeEntities, analyzeSentiment, classifyText.
    - Keep requests modest to respect quotas and billing.
    """

    BASE = "https://language.googleapis.com/v1/documents"

    def __init__(self, api_key: Optional[str] = None, default_language: Optional[str] = None, sleep_between: float = 0.0):
        self.api_key = api_key or os.getenv("GOOGLE_NLP_API_KEY")
        if not self.api_key:
            raise RuntimeError("Google NLP API key missing. Set GOOGLE_NLP_API_KEY or pass api_key.")
        self.default_language = default_language
        self.sleep_between = float(sleep_between or 0.0)

    def _post(self, path: str, payload: Dict[str, Any]) -> Dict[str, Any]:
        url = f"{self.BASE}:{path}?key={self.api_key}"
        resp = requests.post(url, json=payload, timeout=30)
        if resp.status_code != 200:
            try:
                err = resp.json()
            except Exception:
                err = {"error": {"message": resp.text}}
            raise RuntimeError(f"Google NLP {path} error: {resp.status_code} {err}")
        if self.sleep_between:
            time.sleep(self.sleep_between)
        return resp.json()

    def _doc(self, text: str, language: Optional[str] = None) -> Dict[str, Any]:
        doc: Dict[str, Any] = {"type": "PLAIN_TEXT", "content": text}
        lang = language or self.default_language
        if lang:
            doc["language"] = lang
        return doc

    def analyze_entities(self, text: str, language: Optional[str] = None) -> Dict[str, Any]:
        return self._post("analyzeEntities", {"document": self._doc(text, language), "encodingType": "UTF8"})

    def analyze_entity_sentiment(self, text: str, language: Optional[str] = None) -> Dict[str, Any]:
        return self._post("analyzeEntitySentiment", {"document": self._doc(text, language), "encodingType": "UTF8"})

    def analyze_sentiment(self, text: str, language: Optional[str] = None) -> Dict[str, Any]:
        return self._post("analyzeSentiment", {"document": self._doc(text, language), "encodingType": "UTF8"})

    def classify_text(self, text: str, language: Optional[str] = None) -> Dict[str, Any]:
        # classifyText requires sufficient length; caller should guard.
        return self._post("classifyText", {"document": self._doc(text, language)})


def aggregate_insights(
    texts: List[str],
    client: GoogleNLPClient,
    do_entities: bool = True,
    do_entity_sentiment: bool = True,
    do_sentiment: bool = True,
    do_categories: bool = True,
    language: Optional[str] = None,
    max_chars_per_doc: int = 8000,
) -> Dict[str, Any]:
    """Analyze a list of texts and aggregate insights across documents.

    Returns a dict with keys: entities, entity_sentiment, sentiment, categories
    where each value is an aggregated structure ready to convert to DataFrame.
    """

    from collections import defaultdict

    agg_entities = defaultdict(lambda: {"name": None, "type": None, "salience_sum": 0.0, "mentions": 0})
    agg_entity_sent = defaultdict(lambda: {"name": None, "type": None, "sent_sum": 0.0, "magnitude_sum": 0.0, "mentions": 0})
    categories: Dict[str, float] = defaultdict(float)
    sent_scores: List[Tuple[float, float]] = []  # (score, magnitude)

    for txt in texts:
        if not txt:
            continue
        # Trim to keep payloads small and cost reasonable
        t = txt[:max_chars_per_doc]

        if do_entities:
            try:
                e = client.analyze_entities(t, language=language)
                for ent in e.get("entities", []) or []:
                    name = ent.get("name")
                    etype = ent.get("type")
                    sal = float(ent.get("salience", 0.0))
                    key = f"{name}|{etype}"
                    cur = agg_entities[key]
                    cur["name"] = name
                    cur["type"] = etype
                    cur["salience_sum"] += sal
                    cur["mentions"] += max(1, len(ent.get("mentions", []) or [None]))
            except Exception:
                pass

        if do_entity_sentiment:
            try:
                es = client.analyze_entity_sentiment(t, language=language)
                for ent in es.get("entities", []) or []:
                    name = ent.get("name")
                    etype = ent.get("type")
                    sent = float((ent.get("sentiment") or {}).get("score", 0.0))
                    mag = float((ent.get("sentiment") or {}).get("magnitude", 0.0))
                    key = f"{name}|{etype}"
                    cur = agg_entity_sent[key]
                    cur["name"] = name
                    cur["type"] = etype
                    cur["sent_sum"] += sent
                    cur["magnitude_sum"] += mag
                    cur["mentions"] += max(1, len(ent.get("mentions", []) or [None]))
            except Exception:
                pass

        if do_sentiment:
            try:
                s = client.analyze_sentiment(t, language=language)
                doc_s = s.get("documentSentiment") or {}
                score = float(doc_s.get("score", 0.0))
                mag = float(doc_s.get("magnitude", 0.0))
                sent_scores.append((score, mag))
            except Exception:
                pass

        if do_categories:
            # classifyText requires larger docs; combine title+abstract usually suffices
            try:
                c = client.classify_text(t, language=language)
                for cat in c.get("categories", []) or []:
                    name = cat.get("name")
                    conf = float(cat.get("confidence", 0.0))
                    if name:
                        categories[name] += conf
            except Exception:
                pass

    # Prepare results
    ents = list(agg_entities.values())
    ents.sort(key=lambda d: (d["salience_sum"], d["mentions"]), reverse=True)

    entsent = []
    for v in agg_entity_sent.values():
        m = max(1, v["mentions"])
        avg_sent = v["sent_sum"] / m
        avg_mag = v["magnitude_sum"] / m
        entsent.append({**v, "avg_sentiment": avg_sent, "avg_magnitude": avg_mag})
    entsent.sort(key=lambda d: (abs(d["avg_sentiment"]) * d["avg_magnitude"], d["mentions"]), reverse=True)

    cat_list = [{"category": k, "confidence_sum": v} for k, v in categories.items()]
    cat_list.sort(key=lambda d: d["confidence_sum"], reverse=True)

    if sent_scores:
        avg_score = sum(s for s, _ in sent_scores) / len(sent_scores)
        avg_mag = sum(m for _, m in sent_scores) / len(sent_scores)
    else:
        avg_score = 0.0
        avg_mag = 0.0

    return {
        "entities": ents,
        "entity_sentiment": entsent,
        "sentiment": {"avg_score": avg_score, "avg_magnitude": avg_mag, "n": len(sent_scores)},
        "categories": cat_list,
    }

