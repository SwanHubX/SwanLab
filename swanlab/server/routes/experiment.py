from fastapi import APIRouter

router = APIRouter()


@router.get("")
def get():
    return {"test": "test"}
