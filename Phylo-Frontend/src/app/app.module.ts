import { NgModule } from '@angular/core';
import { BrowserModule } from '@angular/platform-browser';
import { AppRoutingModule } from './app-routing.module';
import { AppComponent } from './app.component';

// General imports
import { HomeComponent } from './layout/home/home.component';

// GOF app imports
import { GofHolderComponent } from './gof-holder/gof-holder.component';
import { MapViewComponent } from './gof-holder/map-view/map-view.component';
import { FooterComponent } from './layout/footer/footer.component';
import { NavbarComponent } from './layout/navbar/navbar.component';
import { HeaderComponent } from './layout/header/header.component';
import { ContentComponent } from './layout/content/content.component';


@NgModule({
  declarations: [
    AppComponent,
    GofHolderComponent,
    HomeComponent,
    FooterComponent,
    NavbarComponent,
    HeaderComponent,
    ContentComponent,
  ],
  imports: [
    BrowserModule,
    AppRoutingModule
  ],
  providers: [],
  bootstrap: [AppComponent]
})
export class AppModule { }
