> [!outline]- Table of Contents
> - [[Guidelines]] 
# List of commands

| cmd                     | Description                                                     |
| ----------------------- | --------------------------------------------------------------- |
| `rustup`                |                                                                 |
| `rustup update`         |                                                                 |
| `rustup self uninstall` |                                                                 |
| `rustup doc`            | offline documentation (Rust book)                               |
| `rustc`                 | compiler                                                        |
| `cargo new`             | creates a new project dir w/ git, `cargo.toml`, and subdirs     |
| `cargo init`            | creates a `cargo.toml` file                                     |
| `cargo build`           | builds, producing a binary in `target/debug`                    |
| `cargo build --release` | builds w/ optimizations, producing a binary in `target/release` |
| `cargo run`             | downloads dependencies, builds, and executes                    |
| `cargo check`           | checks whether build succeeds, but produces no binary           |
| `cargo doc --open`      | build documentation for project dependencies                    |

# Syntax
## Variables
```rust
use std::io;
use rand::Rng;
use std::cmp::Ordering;

fn main() {
	let mut bananas = 5;  // mut -> mutable, otherwise immutable
	let mut guess = String::new();  // utf-8 encoded text
	
	io::stdin()
		.read_line(&mut guess)   // first returns Stdin
		.expect("Failed read");  // then Result 'enumeration' of possible state 'variant'
		// Result in {Ok, Err}, expect() crashes in latter case
	
	println!("guess: {guess}");

	let x = 5;
	let y = 10;
	println!("x = {x} and y+2={}", y + 2);
}
```
#### Shadowing
```rust
fn main() {
    let x = 5;

    let x = x + 1;

    {
        let x = x * 2;
        println!("The value of x in the inner scope is: {x}");  // 12
    }

    println!("The value of x is: {x}");  // 6
}
```
#### Data Types

| Integer | unsigned | float   |                                                       | Number literals                                                                                                         |
| ------- | -------- | ------- | ----------------------------------------------------- | ----------------------------------------------------------------------------------------------------------------------- |
| `i8`    | `u8`     |         | Decimal                                               | `98_222`                                                                                                                |
| `i16`   | `u16`    |         | Hex                                                   | `0xff`                                                                                                                  |
| `i32`\* | `u32`    | `f32`   | Octal                                                 | `0o77`                                                                                                                  |
| `i64`   | `u64`    | `f64`\* | Binary                                                | `0b1111_0000`                                                                                                           |
| `i128`  | `u128`   |         | Byte `u8` only                                        | `b'A'`                                                                                                                  |
| `isize` | `usize`  |         | `bool`                                                | `{rust}true; false`                                                                                                     |
|         |          |         | `char`                                                |                                                                                                                         |
|         |          |         | **Tuples**<br>unpacking possible<br>also dot indexing | `{rust}let t: (i32, f64, u8) = (1, 2, 3);`                                                                              |
|         |          |         | **Array**<br>fixed<br>bracket indexing                | `{rust}let a = [1, 2, 3, 4, 5];`<br>`{rust}let a: [i32; 3] = [1, 2, 3];`<br>`{rust}let a = [3; 5];  // [3, 3, 3, 3, 3]` |
|         |          |         |                                                       |                                                                                                                         |
##### Enum
- Algebraic Data Types
	- (possibly recursive) Sum Type of Product Types
	- Sum Type is equivalent to a Disjoint Union

#### Arrays
```rust
use std::io;

let a = [1, 2, 3, 4, 5];

let mut index = String::new();

io::stdin()
	.read_line(&mut index)
	.expect("Failed to read");
	
let index: usize = index
	.trim()
	.parse()
	.expect("Not a number");
	
let element = a[index];

```

### References
```rust
fn main() {
	let mut x: Box<i32> = Box::new(1);
	let a: i32 = *x;         // *x reads the heap value, so a = 1
	*x += 1;                 // *x on the left-side modifies the heap value,
	                         //     so x points to the value 2
	
	let r1: &Box<i32> = &x;  // r1 points to x on the stack
	let b: i32 = **r1;       // two dereferences get us to the heap value
	
	let r2: &i32 = &*x;      // r2 points to the heap value directly
	let c: i32 = *r2;    // so only one dereference is needed to read it
}
```