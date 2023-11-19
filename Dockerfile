# Use uma imagem base
FROM python:3.8

# Configurar variáveis de ambiente
ENV PYTHONUNBUFFERED 1

# Instalar dependências
COPY requirements.txt .
RUN pip install --no-cache-dir -r requirements.txt

# Copiar código
COPY . .

# Rodar aplicação
CMD ["python", "main.py"]

