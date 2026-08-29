# CHARMM-GUI → OpenMM → Vast.ai: Simulação Completa de Membrana

## Objetivo
Simular KKRPKP (PrP 23-28, carga +4) e NNRPNP (controle, carga 0) interagindo
com bicamada POPC explícita. Medir inserção, perturbação de membrana, e afinamento.

## Passo 1: CHARMM-GUI (web, gratuito, ~15 min por sistema)

### 1a. Preparar PDB do peptídeo
- Usar os PDBs gerados: `md_output/KKRPKP_wt_minimized.pdb`
- Ou gerar na própria interface do CHARMM-GUI

### 1b. CHARMM-GUI Membrane Builder
1. Ir para: https://charmm-gui.org/?doc=input/membrane.bilayer
2. **Step 1** — Protein/Peptide:
   - Upload o PDB do peptídeo
   - Orientation: colocar ACIMA da membrana (z = +2 nm do centro)
3. **Step 2** — Lipid composition:
   - Upper leaflet: 64 POPC
   - Lower leaflet: 64 POPC
   - (Total: 128 lipídios, ~6×6 nm patch)
4. **Step 3** — System size:
   - Water thickness: 2.25 nm (acima e abaixo)
   - Ion: 0.15 M KCl
5. **Step 4** — Force field:
   - CHARMM36m (recomendado para peptídeos)
6. **Step 5** — Equilibration:
   - Usar defaults (6-step protocol)
7. **Step 6** — Input generation:
   - Selecionar **OpenMM**
   - Download o `.tgz`

### 1c. Repetir para cada peptídeo
- KKRPKP (wt, charge +4)
- NNRPNP (neutral, charge 0)
- KKRPKN (charge +3, optional)
- KKRPNP (charge +3, different position, optional)

## Passo 2: Preparar scripts para Vast.ai

### 2a. Estrutura de arquivos
```
membrane_sim/
├── system_KKRPKP/       # do CHARMM-GUI
│   ├── step5_input.pdb
│   ├── step5_input.psf   
│   ├── toppar/
│   └── openmm/
│       ├── step6.6_equilibration.py  # do CHARMM-GUI
│       └── step7_production.py       # do CHARMM-GUI
├── system_NNRPNP/       # do CHARMM-GUI
│   └── (mesma estrutura)
├── run_production.py     # nosso script customizado
├── analyze.py            # análise pós-simulação
└── requirements.txt
```

### 2b. No Vast.ai
1. Criar instância: A40 ($0.20-0.35/h)
2. Template: pytorch/pytorch:latest ou nvidia/cuda:12.2.2
3. Upload via `vastai copy` ou rsync
4. Rodar equilibração + produção

## Passo 3: Produção

### Parâmetros
- Timestep: 2 fs
- Ensemble: NPT (constante P e T)
- Temperature: 310 K (fisiológica)
- Pressure: 1 atm, semi-isotropic (xy separado de z)
- Cutoff: 1.2 nm
- PME para eletrostáticos
- SHAKE/LINCS para bonds com H
- Production: 3 × 300 ns por peptídeo
- Save: a cada 10 ps (coordinates) + a cada 1 ps (energias)
- Checkpoint: a cada 10 ns

### Tempo estimado (A40)
- Equilibração (6 steps): ~30 min
- Produção (300 ns): ~5-7 horas
- Total por peptídeo (3 réplicas): ~18 horas
- Total (2 peptídeos × 3 réplicas): ~36 horas
- Custo Vast.ai (~$0.30/h): ~$11

## Passo 4: Análise

### Métricas principais
1. **Profundidade de inserção** — z do COM do peptídeo vs z do centro da membrana
2. **Espessura da membrana** — distância P-P (fósforo upper → lower leaflet)
3. **Ordem lipídica** — SCD dos acils no raio de 1 nm do peptídeo vs controle
4. **SASA** — área acessível ao solvente do peptídeo
5. **Contatos peptídeo-lipídio** — número e tipo (headgroup vs acyl)
6. **Perfil de densidade** — distribuição de massa ao longo de z
7. **Thinning map** — mapa 2D da espessura da membrana (revela perturbação local)

### v5 Predictions (o que esperar)
- KKRPKP: inserção parcial na camada externa (carpet model)
- KKRPKP: afinamento local da membrana (~0.5-1 nm)
- KKRPKP: aumento de contatos com headgroups fosfato
- NNRPNP: permanece na fase aquosa, sem inserção significativa
- A diferença deve ser proporcional à carga
