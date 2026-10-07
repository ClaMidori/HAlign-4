import sys, subprocess, os

input_fasta = sys.argv[1]
output_fasta = sys.argv[2]
chunk_size = 10000 # Fatias de 10kb

# 1. Lê o FASTA gigante (pegando apenas o ID principal)
seqs = {}
order = []
with open(input_fasta, 'r') as f:
    name = ""
    for line in f:
        line = line.strip()
        if line.startswith(">"):
            # Pega apenas a primeira palavra antes do espaço (ex: '>NC_018495.1')
            name = line.split()[0] 
            order.append(name)
            seqs[name] = ""
        else:
            seqs[name] += line

max_len = max(len(s) for s in seqs.values())
aligned_seqs = {name: "" for name in order}

# 2. Fatiamento e Execução
for start in range(0, max_len, chunk_size):
    end = start + chunk_size
    chunk_in = f"tmp_in_{start}.fa"
    chunk_out = f"tmp_out_{start}.fa"

    # Escreve o bloquinho
    with open(chunk_in, 'w') as f:
        for name in order:
            s = seqs[name][start:end]
            if s:
                f.write(f"{name}\n{s}\n")

    # Roda o abPOA apenas no bloquinho
    subprocess.run(["abpoa", chunk_in, "-r", "1"], stdout=open(chunk_out, 'w'))

    # Lê o resultado e costura
    with open(chunk_out, 'r') as f:
        curr_name = ""
        for line in f:
            line = line.strip()
            if line.startswith(">"):
                # Garante que vai pegar só o ID também na saída
                curr_name = line.split()[0]
            else:
                aligned_seqs[curr_name] += line

    # Limpa os arquivos temporários do bloquinho
    os.remove(chunk_in)
    os.remove(chunk_out)

# 3. Salva o alinhamento gigante final
with open(output_fasta, 'w') as f:
    for name in order:
        f.write(f"{name}\n{aligned_seqs[name]}\n")