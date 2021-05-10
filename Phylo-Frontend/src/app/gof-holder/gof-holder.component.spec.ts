import { ComponentFixture, TestBed } from '@angular/core/testing';

import { GofHolderComponent } from './gof-holder.component';

describe('GofHolderComponent', () => {
  let component: GofHolderComponent;
  let fixture: ComponentFixture<GofHolderComponent>;

  beforeEach(async () => {
    await TestBed.configureTestingModule({
      declarations: [ GofHolderComponent ]
    })
    .compileComponents();
  });

  beforeEach(() => {
    fixture = TestBed.createComponent(GofHolderComponent);
    component = fixture.componentInstance;
    fixture.detectChanges();
  });

  it('should create', () => {
    expect(component).toBeTruthy();
  });
});
