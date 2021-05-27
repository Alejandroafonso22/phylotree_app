import { ComponentFixture, TestBed } from '@angular/core/testing';

import { Molecule3DComponent } from './molecule3-d.component';

describe('Molecule3DComponent', () => {
  let component: Molecule3DComponent;
  let fixture: ComponentFixture<Molecule3DComponent>;

  beforeEach(async () => {
    await TestBed.configureTestingModule({
      declarations: [ Molecule3DComponent ]
    })
    .compileComponents();
  });

  beforeEach(() => {
    fixture = TestBed.createComponent(Molecule3DComponent);
    component = fixture.componentInstance;
    fixture.detectChanges();
  });

  it('should create', () => {
    expect(component).toBeTruthy();
  });
});
