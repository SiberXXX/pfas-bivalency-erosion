import urllib.request, json, gzip, time, sys, io, os

sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8')
os.chdir(os.path.dirname(os.path.abspath(__file__)))

# Extract unique gene IDs
genes = []
with gzip.open('../geo_GSE254408/GSE254408_all_samples_raw_counts.txt.gz', 'rt') as f:
    f.readline()  # comment
    f.readline()  # header
    for line in f:
        genes.append(line.split('\t')[0])

print(f'Total genes to map: {len(genes)}')

# Use mygene.info to batch convert (POST with JSON body)
mapping = {}
chunk_size = 1000

for i in range(0, len(genes), chunk_size):
    chunk = genes[i:i+chunk_size]
    url = 'https://mygene.info/v3/gene'
    payload = json.dumps({
        "ids": chunk,
        "scopes": "ensembl.gene",
        "fields": "symbol",
        "species": "human"
    }).encode('utf-8')

    req = urllib.request.Request(url, data=payload, method='POST',
        headers={
            'Content-Type': 'application/json',
            'Accept': 'application/json'
        })

    try:
        resp = urllib.request.urlopen(req, timeout=60)
        results = json.loads(resp.read())
        for r in results:
            if isinstance(r, dict) and 'symbol' in r and 'query' in r:
                mapping[r['query']] = r['symbol']
        batch_num = i // chunk_size + 1
        total_batches = (len(genes) + chunk_size - 1) // chunk_size
        print(f'Batch {batch_num}/{total_batches}: mapped {len(mapping)} so far')
    except Exception as e:
        print(f'Error at batch {i//chunk_size + 1}: {e}')
        # Try smaller chunks
        for j in range(i, min(i + chunk_size, len(genes)), 200):
            sub = genes[j:j+200]
            payload2 = json.dumps({
                "ids": sub,
                "scopes": "ensembl.gene",
                "fields": "symbol",
                "species": "human"
            }).encode('utf-8')
            req2 = urllib.request.Request(url, data=payload2, method='POST',
                headers={'Content-Type': 'application/json', 'Accept': 'application/json'})
            try:
                resp2 = urllib.request.urlopen(req2, timeout=60)
                results2 = json.loads(resp2.read())
                for r in results2:
                    if isinstance(r, dict) and 'symbol' in r and 'query' in r:
                        mapping[r['query']] = r['symbol']
            except Exception as e2:
                print(f'  Sub-batch error: {e2}')
            time.sleep(0.5)
    time.sleep(0.3)

print(f'\nTotal mapped: {len(mapping)} / {len(genes)}')

# Save
with open('../geo_GSE254408/ensembl_to_symbol.tsv', 'w') as f:
    f.write('ensembl_id\tsymbol\n')
    for eid, sym in sorted(mapping.items()):
        f.write(f'{eid}\t{sym}\n')
print('Saved to ../geo_GSE254408/ensembl_to_symbol.tsv')
