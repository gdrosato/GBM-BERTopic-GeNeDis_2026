from pathlib import Path

retrieved_path = Path('data/manifests/retrieved_pmids_2026-08-27.txt')
final_path = Path('data/manifests/final_corpus_pmids.txt')

retrieved = [x.strip() for x in retrieved_path.read_text().splitlines() if x.strip()]
final = [x.strip() for x in final_path.read_text().splitlines() if x.strip()]

assert len(retrieved) == 55446, len(retrieved)
assert len(final) == 54050, len(final)
assert len(set(retrieved)) == len(retrieved), 'Duplicate PMIDs in retrieval manifest'
assert len(set(final)) == len(final), 'Duplicate PMIDs in final-corpus manifest'
assert all(x.isdigit() for x in retrieved), 'Non-numeric PMID in retrieval manifest'
assert all(x.isdigit() for x in final), 'Non-numeric PMID in final-corpus manifest'
assert retrieved == sorted(retrieved, key=int), 'Retrieval manifest is not numerically sorted'
assert final == sorted(final, key=int), 'Final-corpus manifest is not numerically sorted'

retrieved_set = set(retrieved)
final_set = set(final)
assert final_set < retrieved_set, 'Final corpus should be a strict subset of retrieval'
assert len(retrieved_set - final_set) == 1396
assert '27454254' in retrieved_set
assert '27454254' not in final_set

print('Corpus PMID manifest checks: OK')
