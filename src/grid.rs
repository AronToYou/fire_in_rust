use std::ops::{Add, Sub, Mul, Neg, AddAssign};
use std::fmt;

pub trait GridDisp {
    // fn dims(&self) -> (usize, usize);
    fn at(&self, x: usize, y: usize) -> &dyn fmt::Display;
}

impl<T, const NX: usize, const NY: usize> GridDisp for [[T; NY]; NX] where T: fmt::Display {
    // fn dims(&self) -> (usize, usize) { (NX, NY) }
    fn at(&self, x: usize, y: usize) -> &dyn fmt::Display {
        &self[x][y] as &dyn fmt::Display
    }
}

pub trait Linterp: Add<Output = Self> + Mul<f32, Output = Self> + Copy {}
impl<T> Linterp for T where T: Add<Output = Self> + Mul<f32, Output = Self> + Copy {}

#[repr(C)]
#[derive(Clone, Copy, Debug)]
pub struct P<T>(pub T, pub T);

impl<T> Add for P<T> where T: Add<Output = T> + Copy {
    type Output = Self;
    fn add(self, rhs: Self) -> Self {
        P(self.0 + rhs.0, self.1 + rhs.1)
    }
}

impl<T> Sub for P<T> where T: Sub<Output = T> + Copy {
    type Output = Self;
    fn sub(self, rhs: Self) -> Self {
        P(self.0 - rhs.0, self.1 - rhs.1)
    }
}

impl<T> Mul<T> for P<T> where T: Mul<Output = T> + Copy {
    type Output = Self;
    fn mul(self, rhs: T) -> Self {
        P(self.0 * rhs, self.1 * rhs)
    }
}

impl<T> Neg for P<T> where T: Neg<Output = T> + Copy {
    type Output = Self;
    fn neg(self) -> Self {
        P(-self.0, -self.1)
    }
}

impl<T> AddAssign for P<T> where T: Add<Output = T> + Copy {
    fn add_assign(&mut self, other: Self) {
        *self = P(self.0 + other.0, self.1 + other.1)
    }
}

impl From<P<usize>> for P<f32> {
    fn from(p: P<usize>) -> Self {
        P(p.0 as f32, p.1 as f32)
    }
}

impl P<f32> {
    pub fn floor(&self) -> P<usize> {
        P(self.0.floor() as usize, self.1.floor() as usize)
    }
}

impl fmt::Display for P<f32> {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        let p = f.precision().unwrap_or(2);
        write!(f, "({:.p$},{:.p$})", self.0, self.1, p = p)
    }
}

// impl<const N: usize> fmt::Display for &[P<f32>; N] {
//     fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
//         let p = f.precision().unwrap_or(2);
//         for x in 0..N {
//             write!(f, "({:.p$} "), self[x].0, p = p);
//         }
//         for x in 0..N {
//             write!(f, " {:.p$})"), self[x].1, p = p);
//         }
//     }
// }

