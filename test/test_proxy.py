import random
from fastapi.responses import JSONResponse
from fastapi import FastAPI
import uvicorn

app = FastAPI()


@app.get("/api/test")
async def root():
    random_number = random.randint(1, 30)
    return JSONResponse({"data": random_number}, status_code=200)


if __name__ == "__main__":
    uvicorn.run(app, host="0.0.0.0", port=10101)
