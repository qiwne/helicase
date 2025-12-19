pub trait Chunk: Default {}

pub trait Lexer {
    type Input;

    fn input(&self) -> &Self::Input;
}
