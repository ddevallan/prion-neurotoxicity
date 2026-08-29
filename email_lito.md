Assunto: Compilação de pesquisa sobre mecanismo e tratamentos experimentais para DCJ — possíveis direções pouco exploradas

---

Prezados,

Meu nome é Allan, sou engenheiro de software. Não sou médico nem pesquisador, e quero deixar isso claro desde o início. Escrevo porque soube do diagnóstico do Lito e, ao longo desta semana, me dediquei intensamente a entender a doença de Creutzfeldt-Jakob — o mecanismo, o que está sendo testado, o que falhou e por quê.

O que fiz foi uma compilação de pesquisa explorando não apenas a literatura priônica direta, mas também campos adjacentes — biologia de peptídeos antimicrobianos, biofísica de membranas, separação de fases líquido-líquido, vias de morte celular — que descrevem a mesma física mas não se citam entre si. Usei inteligência artificial (Claude, da Anthropic) como parceiro de pesquisa para buscar e cruzar informação de forma paralela e intensiva. Toda a pesquisa, referências e análises estão abertas em:

**https://github.com/ddevallan/prion-neurotoxicity**

Escrevo porque encontrei direções que podem ser relevantes para a equipe médica que está cuidando do Lito. Não estou oferecendo diagnóstico nem tratamento — estou compartilhando informação compilada que pode valer a pena ser avaliada por especialistas.

---

## O que encontrei que pode ser imediatamente relevante

### 1. Memantina e antagonistas de NMDA — história longa, nunca adequadamente perseguida

A memantina (Namenda/Ebix) é um antagonista de NMDA aprovado pela FDA para Alzheimer. Existe uma trilha de evidências sobre antagonistas de NMDA na doença priônica que se estende por mais de 50 anos, mas nunca foi adequadamente seguida:

**a) Amantadina em pacientes com DCJ (1971–1984)**
A amantadina (antagonista fraco de NMDA, ~10× menos potente que memantina) foi testada em pacientes com DCJ nos anos 1970. Sanders & Dunn 1973 relataram dois casos: um paciente teve melhora notável por 2 meses; o outro **"aparentou estar curado" no follow-up de 30 meses.** Uma revisão sistemática posterior (Stewart 2008) encontrou **4 de 8 relatos com benefício.** Na época, não se sabia que amantadina era antagonista de NMDA — foi testada como antiviral. A conexão NMDA nunca foi feita retroativamente.

**b) Memantina protege in vitro (1993)**
Müller et al. (European Journal of Pharmacology 1993) demonstraram que memantina, um análogo, e MK-801 preveniram a morte neuronal induzida por PrPSc em culturas corticais de rato.

**c) Memantina testada in vivo — resultado modesto mas significativo (2008)**
Riemer et al. 2008 (Alemanha) testaram memantina em camundongos infectados com scrapie 139A: 30 mg/kg/dia oral, iniciada no dia 100 pós-infecção (tarde). **Resultado: +8% de extensão de sobrevida (196±11.7 vs 181±7.2 dias, p<0.01).** Foi o melhor resultado entre 5 drogas testadas. Porém: dose suprafisiológica (6-24× a dose humana), início tardio, n=8, nunca seguido com timing otimizado.

**d) Memantina resgata modelo genético (2025)**
JCI 2025 (doi:10.1172/JCI186432): camundongos Prnp G92N com influxo persistente de Ca²⁺ via NMDA, beading dendrítico e convulsões — tudo resgatado por memantina.

**e) PrPC normalmente inibe NMDA (2008)**
Khosravani et al. 2008 mostraram que PrPC inibe subunidades NR2D do NMDA. Quando PrPC é sequestrada pela conversão priônica, essa inibição é perdida — abrindo um segundo mecanismo de excitotoxicidade.

**f) Não é só NMDA — receptores AMPA também estão envolvidos**
Bhatt et al. 2020 (PLoS Pathogens) mostraram que PrP mutante retém subunidades GluA2, deixando receptores AMPA permeáveis a cálcio na sinapse (10× mais morte neuronal). **Perampanel** (Fycompa), um antagonista de AMPA aprovado para epilepsia, poderia cobrir esse eixo.

**O quadro completo:** a toxicidade priônica envolve excitotoxicidade por DOIS eixos (NMDA + AMPA). Existem drogas aprovadas para ambos. A amantadina (fraca) mostrou benefício em pacientes nos anos 1970. A memantina (forte) mostrou +8% em camundongo com design subótimo. Nenhuma foi testada com timing otimizado ou em combinação.

**Perfil de segurança:** memantina penetra o LCR (52% do plasma), biodisponibilidade oral de 100%, meia-vida 60-80h, sem interações com trazodona, lítio ou ASO intratecal. Sem contraindicação que impeça uso em paciente com DCJ.

**Referências:**
- Sanders WL, Dunn TL. 1973 (PMC494412) — amantadina em DCJ
- Müller WE et al. Eur J Pharmacol 1993 (PMID 7901042) — memantina in vitro
- Riemer et al. 2008 — memantina in vivo, +8%
- JCI 2025, doi:10.1172/JCI186432 — resgate do G92N
- Khosravani et al. 2008 (PMC2364707) — PrPC inibe NR2D
- Bhatt et al. PLoS Pathog 2020 — AMPA, GluA2
- Resenberger et al. 2011 (PMC3098494) — mecanismo NMDA compartilhado

**Ressalva:** memantina não para a replicação priônica. O +8% de Riemer 2008 é consistente com o esperado para intervenções downstream (Minikel documenta 4-19%). Não é cura — mas pode proteger neurônios enquanto intervenções que reduzem PrP (ASO/siRNA) fazem efeito.

---

### 2. Ensaios clínicos ativos — dois trials recrutando agora

Existem dois ensaios clínicos ativos para doença priônica em agosto de 2026:

**a) PrProfile (ION717, Ionis Pharmaceuticals)**
- Oligonucleotídeo antisense (ASO) que reduz a produção de PrPC
- Fase 1/2a, 56 pacientes sintomáticos inscritos
- Terceiro esquema de dose adicionado em março de 2026
- Administração intratecal (punção lombar)
- Conclusão primária estimada para fevereiro de 2027

**b) PRiSM (Broad Institute / UMass Chan / Regeneron)**
- siRNA divalente, dose única intracerebroventricular
- Fase 1, 15 pacientes sintomáticos
- IND do FDA aprovado em março de 2025, recrutando desde abril de 2026
- Dados pré-clínicos: 49% de redução de PrP cerebral, 64% de extensão de sobrevida em camundongos mesmo quando iniciado após início dos sintomas
- PI: Eric Minikel (pesquisador que perdeu a esposa para a doença e se tornou cientista para buscar cura)

A equipe médica provavelmente já conhece esses trials. Se não, vale verificar elegibilidade.

---

### 3. Outras drogas já aprovadas que podem ter relevância

**a) Lítio (dose baixa)**
- Lítio em dose baixa melhorou sobrevida, reduziu vacuolização e perda neuronal em camundongos infectados com príon
- Mecanismo: indução de autofagia via via independente de mTOR
- Já aprovado e disponível
- Referência: PubMed 30135493

**b) Trazodona**
- Antidepressivo que inibe a via PERK/eIF2α (resposta a proteínas mal-dobradas)
- Preveniu neurodegeneração em camundongos infectados com príon quando administrada oralmente
- Restaurou tradução proteica e protegeu sinapses e mitocôndrias
- Referência: Halliday et al. Brain 2017

**Ressalva sobre todas essas drogas:** nenhuma delas para a replicação priônica. Elas protegem neurônios contra as consequências da doença. A única classe que demonstrou alterar fundamentalmente o curso é a redução de PrP (os ASOs/siRNAs dos trials acima).

---

### 4. O que NÃO funciona (para evitar buscas em direções erradas)

A pesquisa que compilei também mapeou tratamentos que falharam e por quê:

- **Quinacrina:** zero benefício em trial clínico; gerou príons resistentes à droga
- **Doxiciclina:** nunca teve dados pré-clínicos positivos em modelos relevantes
- **Pentosana polissulfato (PPS):** requer infusão intraventricular com complicações; sem benefício claro em humanos
- **Anticorpo PRN100 (UCL):** 6 pacientes tratados, zero sinal de eficácia

O padrão transversal das falhas é: **nada funciona depois do início clínico avançado, exceto (parcialmente) a redução de PrP.** E atacar a PrPSc diretamente é atacar algo que pode evoluir resistência; atacar a PrPC (o substrato do hospedeiro) é atacar algo que não pode escapar.

---

### 5. O modelo que emergiu da pesquisa (para os médicos/pesquisadores)

O modelo completo está no repositório (docs/MODEL_v5.md), mas em resumo:

A doença priônica pode ser entendida como "fogo amigo" — a proteína priônica normal (PrPC) tem uma região N-terminal que funciona como peptídeo antimicrobiano (confirmado por Nature Microbiology, junho 2026 — "prionins"). Durante a conversão priônica, essa região é liberada na membrana neuronal e faz ao neurônio o que evoluiu para fazer a bactérias: perturba a membrana.

A especificidade da doença priônica não está no mecanismo de toxicidade (que é genérico — qualquer peptídeo catiônico suficiente faz o mesmo), mas na velocidade de replicação: príons se replicam ~1000× mais rápido que tau ou alfa-sinucleína, o que explica por que a doença mata em meses enquanto Alzheimer leva décadas.

O protetor natural mais poderoso da espécie — a variante G127V, dos Fore de Papua Nova Guiné que sobreviveram ao kuru — funciona mantendo os condensados proteicos em estado líquido, impedindo a transição para fibrilas sólidas. Ele não toca no mecanismo de toxicidade; bloqueia a conversão.

---

### 6. O que a combinação ideal seria (segundo o modelo)

Nenhum tratamento isolado é suficiente depois do início dos sintomas. O modelo sugere uma combinação:

1. **Redução de PrPC** (ASO ou siRNA — os trials acima) — corta a fonte
2. **Memantina** (já aprovada) — protege neurônios na janela enquanto o ASO faz efeito
3. **Trazodona** (já aprovada) — protege a maquinaria de tradução
4. **Lítio baixa dose** (já aprovado) — induz autofagia, limpa agregados

Os itens 2, 3 e 4 são drogas já aprovadas e disponíveis. O item 1 é experimental e requer inscrição em trial.

**Ninguém testou essa combinação.** Cada trial testa uma droga isolada. Mas sob o modelo, a combinação não é luxo — é necessidade, porque cada componente age num passo diferente da cascata.

---

## Limitações importantes

- **Não sou especialista.** Tudo que encontrei precisa ser avaliado por profissionais qualificados.
- **Nada aqui é revisado por pares.** É uma investigação exploratória, não uma revisão sistemática.
- **As referências foram verificadas contra resumos e trechos**, mas detalhes bibliográficos devem ser confirmados.
- **O modelo pode estar errado.** Ele ficou elegante demais, o que é suspeito. As limitações estão documentadas honestamente no repositório.
- **Nada aqui substitui orientação médica.**

---

## Contato e repositório

Todo o material — modelo, referências (~150 artigos), análises computacionais, mapa terapêutico, experimentos propostos — está aberto em:

**https://github.com/ddevallan/prion-neurotoxicity**

Se houver qualquer pesquisador ou médico na equipe que queira discutir os achados, conversar sobre as referências, ou apontar onde estou errado, ficarei genuinamente grato. Meu objetivo não é estar certo — é que a informação chegue a quem pode avaliá-la.

Allan
allanmfx@gmail.com
