#include <iostream>
#include <iomanip>
#include <cmath>
#include <fstream>
#include <string>
#include <cstdio>

using namespace std;

const double PI = 3.14159265358979323846;
const double VMIN = 500.0;
const double DELTA = 1e-5;
const int MAX_ITER = 100000;
const double TOL_GRAD = 1e-4;
const double TOL_MOV = 1e-9;
const double ALPHA_MIN = 1e-16;

// Custo: base + parede lateral + teto
double custo(double r, double h) {
    return 350.0 * PI * r * r + 200.0 * PI * r * h;
}

double volume(double r, double h) {
    return PI * r * r * h;
}

// g(r,h) <= 0 representa a restricao de volume
double restricao(double r, double h) {
    return VMIN - volume(r, h);
}

double penalizada(double r, double h, double mu) {
    double g = restricao(r, h);
    if (g < 0.0) g = 0.0;
    return custo(r, h) + mu * g * g;
}

// Diferencas finitas progressivas para a funcao penalizada geral.
double derivadaR(double r, double h, double mu) {
    return (penalizada(r + DELTA, h, mu) - penalizada(r, h, mu)) / DELTA;
}

double derivadaH(double r, double h, double mu) {
    return (penalizada(r, h + DELTA, mu) - penalizada(r, h, mu)) / DELTA;
}

// A solucao de custo minimo tem volume ativo: V = VMIN.
// Portanto, h pode ser calculada a partir de r, mantendo a viabilidade:
// h(r) = VMIN / (PI*r^2). A funcao abaixo e o custo penalizado nessa fronteira.
double penalizadaViavel(double r, double mu) {
    if (r <= 0.0) return 1e100;
    double h = VMIN / (PI * r * r);
    return penalizada(r, h, mu);
}

// Derivada numerica da funcao reduzida, mantendo V = VMIN.
double derivadaReduzida(double r, double mu) {
    return (penalizadaViavel(r + DELTA, mu) -
            penalizadaViavel(r, mu)) / DELTA;
}

struct Resultado {
    double raio;
    double altura;
    double volumeFinal;
    double custoFinal;
    double normaGradiente;
    int iteracoes;
    bool convergiu;
    string motivo;
};

Resultado otimizar(double mu, double alphaInicial,
                   const string& arquivoHistorico) {
    Resultado res;
    double r = 3.5;
    double h = VMIN / (PI * r * r);
    double alpha = alphaInicial;
    double pAtual, gr, norma;
    double nr, nh, np, alphaTeste;
    double deslocamento;
    int i, tentativa;
    bool achou;

    res.convergiu = false;
    res.iteracoes = 0;
    res.normaGradiente = 0.0;
    res.motivo = "Limite maximo de iteracoes atingido";

    ofstream hist(arquivoHistorico.c_str());
    if (!hist.is_open()) {
        res.raio = r;
        res.altura = h;
        res.volumeFinal = volume(r, h);
        res.custoFinal = custo(r, h);
        res.motivo = "Erro ao criar arquivo historico CSV";
        return res;
    }

    hist << "Iteracao,Raio,Altura,Volume,Custo,Penalizada,"
         << "DerivadaReduzida,NormaGradiente,Alpha\n";

    for (i = 0; i < MAX_ITER; ++i) {
        h = VMIN / (PI * r * r);
        pAtual = penalizada(r, h, mu);
        gr = derivadaReduzida(r, mu);
        // Norma do gradiente reduzido: a direcao respeita a restricao V=500.
        norma = fabs(gr);
        res.normaGradiente = norma;

        hist << setprecision(12)
             << i << "," << r << "," << h << ","
             << volume(r, h) << "," << custo(r, h) << ","
             << pAtual << "," << gr << "," << norma << "," << alpha << "\n";

        if (norma < TOL_GRAD) {
            res.convergiu = true;
            res.motivo = "Convergencia: derivada reduzida abaixo da tolerancia";
            break;
        }

        // Busca em linha na variavel r. A altura e recalculada a cada teste,
        // garantindo que todo candidato tenha volume igual a 500 m^3.
        alphaTeste = alpha;
        achou = false;
        for (tentativa = 0; tentativa < 100; ++tentativa) {
            nr = r - alphaTeste * gr;
            if (nr > 0.0) {
                nh = VMIN / (PI * nr * nr);
                np = penalizada(nr, nh, mu);
                if (np < pAtual) {
                    achou = true;
                    break;
                }
            }
            alphaTeste *= 0.5;
            if (alphaTeste < ALPHA_MIN) break;
        }

        if (!achou) {
            // Para perto do minimo, reduza o passo e tente novamente.
            alpha *= 0.1;
            if (alpha < 1e-18) {
                res.motivo = "Parada: passo muito pequeno sem descida numerica";
                break;
            }
            continue;
        }

        deslocamento = fabs(nr - r);
        r = nr;
        alpha = alphaTeste;
        res.iteracoes = i + 1;

        if (deslocamento < TOL_MOV && norma < 10.0 * TOL_GRAD) {
            res.convergiu = true;
            res.motivo = "Convergencia: deslocamento e derivada reduzida pequenos";
            break;
        }
    }

    if (i >= MAX_ITER) res.iteracoes = MAX_ITER;
    hist.close();

    res.raio = r;
    res.altura = VMIN / (PI * r * r);
    res.volumeFinal = volume(res.raio, res.altura);
    res.custoFinal = custo(res.raio, res.altura);
    if (res.convergiu && fabs(res.volumeFinal - VMIN) > 1e-6) {
        res.convergiu = false;
        res.motivo = "Parada: erro numerico na restricao de volume";
    }
    return res;
}

void mostrarResultado(double mu, double alpha, const Resultado& res) {
    cout << "\n----------------------------------------\n";
    cout << "RESULTADO DA OTIMIZACAO\n";
    cout << "----------------------------------------\n";
    cout << setprecision(10);
    cout << "Mu: " << mu << "\n";
    cout << "Alpha inicial: " << alpha << "\n";
    cout << "Raio final: " << res.raio << " m\n";
    cout << "Altura final: " << res.altura << " m\n";
    cout << "Volume final: " << res.volumeFinal << " m^3\n";
    cout << "Custo final: R$ " << res.custoFinal << "\n";
    cout << "Norma do gradiente: " << res.normaGradiente << "\n";
    cout << "Iteracoes: " << res.iteracoes << "\n";
    cout << "Convergencia do metodo: "
         << (res.convergiu ? "SIM" : "NAO") << "\n";
    cout << "Motivo de parada: " << res.motivo << "\n";
    cout << "Restricao de volume: "
         << (res.volumeFinal >= VMIN ? "ATENDIDA" : "NAO ATENDIDA")
         << "\n";
    cout << "----------------------------------------\n";
}

void testeIndividual() {
    double mu, alpha;
    char nome[160];

    cout << "\nDigite mu (ex.: 10000): ";
    cin >> mu;
    cout << "Digite alpha inicial (ex.: 0.001): ";
    cin >> alpha;

    if (!cin || mu <= 0.0 || alpha <= 0.0) {
        cout << "Erro: digite valores numericos positivos.\n";
        return;
    }

    sprintf(nome, "historico_mu_%.0f_alpha_%.6g.csv", mu, alpha);
    Resultado res = otimizar(mu, alpha, nome);
    mostrarResultado(mu, alpha, res);
    cout << "Historico salvo em: " << nome << "\n";
}

void compararParametros() {
    const double mus[3] = {1000.0, 10000.0, 100000.0};
    const double alphas[2] = {1e-4, 1e-5};
    ofstream resumo("resumo_testes.csv");
    int i, j;

    if (!resumo.is_open()) {
        cout << "Erro ao criar resumo_testes.csv\n";
        return;
    }

    resumo << "Mu,AlphaInicial,Raio,Altura,Volume,Custo,NormaGradiente,"
           << "Iteracoes,Convergiu,Motivo\n";

    for (i = 0; i < 3; ++i) {
        for (j = 0; j < 2; ++j) {
            char nome[160];
            double mu = mus[i];
            double alpha = alphas[j];
            sprintf(nome, "historico_mu_%.0f_alpha_%.6g.csv", mu, alpha);

            cout << "\nExecutando mu = " << mu
                 << ", alpha = " << alpha << "\n";
            Resultado res = otimizar(mu, alpha, nome);
            mostrarResultado(mu, alpha, res);

            resumo << setprecision(12)
                   << mu << "," << alpha << ","
                   << res.raio << "," << res.altura << ","
                   << res.volumeFinal << "," << res.custoFinal << ","
                   << res.normaGradiente << "," << res.iteracoes << ","
                   << (res.convergiu ? "SIM" : "NAO") << ","
                   << res.motivo << "\n";
        }
    }

    resumo.close();
    cout << "\nResumo salvo em resumo_testes.csv\n";
}

int main() {
    int opcao;

    cout << "========================================\n";
    cout << " OTIMIZACAO DO CUSTO DE UM TANQUE\n";
    cout << "========================================\n";
    cout << "Volume minimo: 500 m^3\n\n";
    cout << "1 - Teste individual\n";
    cout << "2 - Comparar diferentes valores de mu e alpha\n";
    cout << "Escolha uma opcao: ";
    cin >> opcao;

    if (opcao == 1) {
        testeIndividual();
    } else if (opcao == 2) {
        compararParametros();
    } else {
        cout << "Opcao invalida.\n";
    }

    cout << "\nPrograma finalizado. Pressione ENTER para sair...";
    cin.ignore(1000, '\n');
    cin.get();
    return 0;
}

