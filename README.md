# SARS-CoV-2 Analysis Toolkit

Este repositorio recopila utilidades utilizadas durante el esfuerzo de secuenciación
SARS-CoV-2 para el Sistema Sanitario Público de Andalucía. Incluye un script para
orquestar la limpieza y agregación de metadatos, una canalización Nextstrain y un
notebook histórico que sirvió como referencia inicial.

## Contenido del repositorio

- `parse_sars-cov2_metadata.py`: versión en script del flujo de trabajo de
  metadatos empleado originalmente en `tercera_ola.ipynb`. Recolecta información
  desde distintas carpetas de red, normaliza identificadores, enriquece con
  datos epidemiológicos y genera salidas listas para Nextstrain y reportes
  web.
- `nextstrain_pipeline.sh`: secuencia de comandos basada en `augur` para filtrar,
  alinear y exportar conjuntos de datos hacia `auspice`.
- `tercera_ola.ipynb`: cuaderno original desde el que se derivó el script. Se
  conserva como documentación del flujo exploratorio inicial.

## Requisitos

### Python

El script de metadatos requiere Python 3.9 o superior y las siguientes
dependencias:

- pandas
- numpy
- openpyxl (para generar los ficheros Excel de control)

Se recomienda crear un entorno virtual y luego instalar las dependencias:

```bash
python -m venv .venv
source .venv/bin/activate
pip install pandas numpy openpyxl
```

### Nextstrain

Para ejecutar la canalización `nextstrain_pipeline.sh` se necesita tener instalados
los componentes del ecosistema Nextstrain (`augur`, `auspice` y las bases de datos
de referencia). Consulte la [documentación oficial de Nextstrain](https://docs.nextstrain.org/projects/augur/en/stable/) para
preparar el entorno.

## Estructura de datos esperada

Los scripts asumen un árbol de directorios que reside, por defecto, en la ruta
`/mnt/NAS/projects/virus_SSPA/COVID19_SSPA/`. Se puede sobrescribir la ubicación
principal definiendo la variable de entorno `COVID19_SSPA_BASE` antes de lanzar
el script:

```bash
export COVID19_SSPA_BASE=/ruta/a/COVID19_SSPA
```

Dentro de este árbol se esperan las subcarpetas `nextstrain_all/data`,
`other_data/1wave_data` y `secuenciacion_salud_publica/metadata` con la misma
organización utilizada durante el proyecto original.

## Ejecución del procesamiento de metadatos

1. Configure el entorno Python y verifique que las rutas de datos sean accesibles.
2. (Opcional) Ajuste la variable `COVID19_SSPA_BASE` si los datos residen en una
   ubicación diferente.
3. Ejecute el script:

   ```bash
   python parse_sars-cov2_metadata.py
   ```

El proceso generará:

- Tablas consolidadas en `analysis_MN908947/other/` dentro de la carpeta de
  secuenciación.
- Ficheros `TSV` agregados para uso web en `analysis_MN908947/reports/web/`.
- Los ficheros `data/auspice_metadata.tsv` y `data/sequence.fasta` que sirven de
  entrada al pipeline de Nextstrain.

Si se detectan inconsistencias (identificadores duplicados, fechas fuera de
rango, etc.), el script imprimirá mensajes de advertencia y dejará ficheros de
apoyo en el directorio de trabajo (por ejemplo, `error_date_lineage.xlsx`).

## Ejecución de la canalización Nextstrain

Una vez generados los ficheros de `data/`, cree las carpetas `results/` y
`auspice/` si no existen. Después ejecute:

```bash
bash nextstrain_pipeline.sh
```

El script aplicará filtrados, reconstruirá árboles filogenéticos, generará
mutaciones ancestrales y exportará un JSON compatible con Auspice. El resultado
final se enlaza como `auspice/SARS-COV-2-latest.json` para facilitar su
visualización.

## Notebook de referencia

El cuaderno `tercera_ola.ipynb` se mantiene como documentación. No forma parte
del flujo automatizado pero es útil para comprender el origen de los procesos
incluidos en el script.

## Contribución

Para extender o adaptar estas utilidades:

1. Abra un nuevo branch en Git.
2. Aplique los cambios y añada pruebas manuales o automáticas según corresponda.
3. Documente cualquier nueva dependencia o ruta esperada en este README.

Las contribuciones que simplifiquen el mantenimiento y reduzcan dependencias
externas son especialmente bienvenidas.
