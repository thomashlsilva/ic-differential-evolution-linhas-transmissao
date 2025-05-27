# ⚡ Estudo dos Campos Elétricos em Linhas de Transmissão Aéreas

Este repositório contém os códigos e modelos desenvolvidos durante a iniciação científica intitulada:

**"Estudo das configurações geométricas das linhas de transmissão aéreas com vistas a obtenção de melhorias no campo elétrico"**

## 📑 Resumo

Este trabalho visa estudar as configurações geométricas de linhas de transmissão aéreas com elevada capacidade de transmissão empregadas no Brasil. O estudo é realizado por meio do levantamento das características físicas e elétricas destes equipamentos.

Foram analisados:

- Efeitos **Joule** e **Corona**.
- Componentes físicos das linhas, como **condutores, isoladores, cabos para-raios** e estruturas de sustentação.

A metodologia inclui:

- Modelagem dos campos elétricos por meio do **Método das Imagens** e do **Método das Imagens Sucessivas**.
- Cálculo do campo elétrico ao **nível do solo** e do **campo elétrico superficial**.
- Formulação de um problema de **otimização restrita**, utilizando o **Algoritmo de Evolução Diferencial (Differential Evolution)** com modificação no parâmetro de amplificação diferencial (*F*).

A validação foi realizada comparando os resultados dos algoritmos com dados da literatura, considerando diferentes configurações (1, 2, 3 e 4 cabos por fase), respeitando as normas técnicas da **NBR 5422** para linhas aéreas de transmissão.

Por fim, desenvolveu-se uma ferramenta computacional capaz de calcular e minimizar os níveis de campo elétrico superficial e ao nível do solo, garantindo conformidade normativa e melhor desempenho elétrico.

---

## 🗂️ Organização do Repositório

```
.
├── misc/           # Scripts auxiliares e arquivos de anotação
├── otimizacao/     # Funções de cálculo dos campos elétricos e do algoritmo Differential Evolution com F modificado
├── solo/           # Modelos de cálculo do campo elétrico ao nível do solo (4 casos tratados)
├── superficial/    # Modelos de cálculo do campo elétrico superficial (5 casos tratados)
├── variacao/       # Scripts para gerar variações nos parâmetros geométricos (uso didático e de exemplificação)
└── README.md       # Este arquivo
```

---

## 🚀 Como Executar

- Os scripts foram desenvolvidos em ambiente **MATLAB/Octave**.
- ✅ **Compatível com GNU Octave**, uma alternativa **livre e de código aberto** ao MATLAB.
- Basta executar os arquivos `.m` localizados nas pastas conforme o objetivo desejado.

### 📌 Exemplos:
- Cálculo de campo elétrico ao nível do solo:
  - Navegue até `solo/` e execute o script referente ao caso de interesse.

- Cálculo do campo elétrico superficial:
  - Acesse `superficial/` e execute o script do cenário desejado.

- Execução do algoritmo de **Evolução Diferencial**:
  - Acesse `otimizacao/` e execute o script principal de otimização, juntamente com as funções auxiliares presentes na mesma pasta.

---

## 📚 Como Citar

```bibtex
@misc{silva2021,
  author  = {Thomás Henrique Lopes Silva},
  title   = {Estudo das configurações geométricas das linhas de transmissão aéreas com vistas a obtenção de melhorias no campo elétrico},
  school  = {Centro Federal de Educação Tecnológica de Minas Gerais},
  year    = {2021},
  address = {Divinópolis, Brasil},
  note    = {Código e dados disponíveis em: \url{https://github.com/thomashlsilva/ic-differential-evolution-linhas-transmissao}}
}
```

---

**Contato:**
✉️ [thomashlsilva@gmail.com](mailto:thomashlsilva@gmail.com)
💻 [@thomashlsilva](https://github.com/thomashlsilva)

---

## 📜 Licença

* Códigos: [![MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)
* Dados/Resultados: [![CC BY 4.0](https://img.shields.io/badge/License-CC_BY_4.0-lightgrey.svg)](LICENSE-CC-BY-4.0)
