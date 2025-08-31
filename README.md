# CancerGeneSignatures
Generador de firmas genéticas asociadas a tipos de cáncer mediante análisis de expresión génica, integración bibliográfica y construcción de redes de interacción. Facilita la identificación de genes clave y su aplicación en investigación oncológica, diagnóstico o pronóstico.

## 🚀 Características Principales

- **Análisis de Genes Normalizadores**: Identificación robusta de genes de referencia para qPCR
- **Validación Estadística Completa**: Múltiples estrategias de validación de genes normalizadores
- **Análisis de Expresión Diferencial**: Detección de genes con expresión alterada
- **Integración Bibliográfica**: Análisis de literatura científica
- **Construcción de Redes**: Análisis de interacciones génicas

## 📋 Módulos Principales

### 🔬 Validación de Genes Normalizadores (`src/validation.py`)
Módulo completo para validar la robustez de genes normalizadores seleccionados:

- **Validación Cruzada K-fold**: Evalúa estabilidad en diferentes subconjuntos
- **Análisis de Bootstrap**: Determina robustez mediante remuestreo
- **Validación de Estabilidad Temporal**: Analiza consistencia a lo largo del tiempo
- **Análisis de Sensibilidad**: Evalúa robustez ante perturbaciones
- **Validación de Correlación Técnica**: Verifica independencia de variables técnicas
- **Análisis de Reproducibilidad**: Valida consistencia inter-experimental

### 📊 Análisis de Datos (`src/data_processing.py`)
Procesamiento y análisis de datos de expresión génica.

### 💬 Mensajería (`src/messaging.py`)
Sistema de notificaciones y reportes.

### 💾 Guardado de Resultados (`src/save_results.py`)
Gestión de resultados y exportación de datos.

## 🛠️ Instalación

1. Crea un entorno virtual (opcional pero recomendado):
   ```bash
   python -m venv .venv
   source .venv/bin/activate
   ```

2. Instala las dependencias necesarias:
   ```bash
   pip install -r requirements.txt
   ```

## 📖 Uso

### Validación de Genes Normalizadores

```python
from src.validation import validacion_completa_genes_normalizadores
import pandas as pd

# Cargar datos
controles_df = pd.read_csv('datos_controles.csv')
muestras_df = pd.read_csv('datos_muestras.csv')

# Genes de referencia a validar
genes_referencia = ['GAPDH', 'ACTB', '18S']

# Ejecutar validación completa
validador = validacion_completa_genes_normalizadores(
    df_controles=controles_df,
    df_muestras=muestras_df,
    genes_referencia=genes_referencia,
    generar_graficos=True,
    guardar_reporte=True
)
```

### Script de Ejemplo

```bash
# Ejecutar ejemplo de validación
python src/ejemplo_validacion.py
```

### Notebook de Validación

```bash
# Abrir notebook de validación
jupyter notebook notebooks/validacion_genes_normalizadores.ipynb
```

## 📊 Estrategias de Validación Implementadas

### 1. Validación Cruzada K-fold
- **Objetivo**: Evaluar estabilidad de genes de referencia
- **Métrica**: Correlación entre estabilidades en train/test
- **Criterio**: Score > 0.7 indica buena estabilidad

### 2. Análisis de Bootstrap
- **Objetivo**: Determinar robustez mediante remuestreo
- **Métrica**: Frecuencia de aparición como gen más estable
- **Criterio**: Frecuencia > 50% indica alta robustez

### 3. Análisis de Sensibilidad
- **Objetivo**: Evaluar robustez ante perturbaciones
- **Métrica**: Cambio porcentual en score de separación
- **Criterio**: Cambio < 10% indica alta robustez

### 4. Validación de Estabilidad Temporal
- **Objetivo**: Analizar consistencia a lo largo del tiempo
- **Métrica**: Correlación con tiempo
- **Criterio**: |Correlación| < 0.3 indica estabilidad temporal

### 5. Validación de Correlación Técnica
- **Objetivo**: Verificar independencia de variables técnicas
- **Métrica**: Correlación con variables técnicas
- **Criterio**: |Correlación| < 0.5 indica independencia

### 6. Análisis de Reproducibilidad
- **Objetivo**: Validad consistencia inter-experimental
- **Métrica**: ICC (Intraclass Correlation Coefficient)
- **Criterio**: ICC > 0.75 indica excelente reproducibilidad

## 📈 Interpretación de Resultados

- **Score Global > 0.8**: Excelente confiabilidad
- **Score Global 0.6-0.8**: Buena confiabilidad  
- **Score Global 0.4-0.6**: Aceptable, requiere monitoreo
- **Score Global < 0.4**: Problemático, reevaluar genes

## 📁 Estructura del Proyecto

```
CancerGeneSignatures/
├── src/
│   ├── validation.py              # Módulo de validación
│   ├── data_processing.py         # Procesamiento de datos
│   ├── messaging.py               # Sistema de mensajería
│   ├── save_results.py            # Guardado de resultados
│   └── ejemplo_validacion.py      # Script de ejemplo
├── notebooks/
│   ├── genes_normalizadores.ipynb # Análisis de genes normalizadores
│   └── validacion_genes_normalizadores.ipynb # Validación completa
├── raw_data/                      # Datos crudos
├── gen-sets_GSEA_MSigDB/          # Conjuntos de genes
└── requirements.txt               # Dependencias
```

## 🔧 Dependencias

- pandas >= 1.3.0
- numpy >= 1.21.0
- matplotlib >= 3.4.0
- seaborn >= 0.11.0
- scipy >= 1.7.0
- scikit-learn >= 1.0.0
- plotly >= 5.0.0
- jupyter >= 1.0.0

## 📚 Referencias

- Vandesompele et al. (2002) - geNorm
- Andersen et al. (2004) - NormFinder  
- Pfaffl et al. (2004) - BestKeeper
- Livak & Schmittgen (2001) - Método ΔΔCt
- Bustin et al. (2009) - Guidelines MIQE

## 🤝 Contribuciones

Las contribuciones son bienvenidas. Por favor, abre un issue o pull request para sugerir mejoras.

## 📄 Licencia

Este proyecto está bajo la Licencia MIT. Ver el archivo `LICENSE` para más detalles.
